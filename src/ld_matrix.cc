#include "ld_matrix.h"

#include <string>
#include <valarray>

#include "bgmg_log.h"
#include "bgmg_parse.h"
#include "plink_ld.h"
#include "ld_matrix_csr.h"

#define LD_MATRIX_FORMAT_VERSION 1

class PosixFile {
public:
  PosixFile(std::string filename, std::string modes) {
    file = fopen(filename.c_str(), modes.c_str());
    if (file == nullptr) BGMG_THROW_EXCEPTION(::std::runtime_error(std::string("Unable to open ") + filename));
  }
  ~PosixFile() { fclose(file); }
  FILE* handle() { return file; }

private:
  FILE* file;
};

void generate_ld_matrix_from_bed_file(std::string bfile, float r2_min, float ldscore_r2min, std::string outfile) {
  LOG << ">generate_ld_matrix_from_bed_file(bfile=" << bfile << ", r2_min=" << r2_min << ", ldscore_r2min=" << ldscore_r2min << ");";
  SimpleTimer timer(-1);

  FamFile fam_file(bfile + ".fam");
  BimFile bim_file(bfile + ".bim");
  bim_file.find_snp_to_index_map();
  const int num_snps = bim_file.size();
  const int num_subj = fam_file.size();

  const int block_size = std::min(8*1024, num_snps);  // handle blocks of up to 8K SNPs
  const int block_elems = block_size * block_size;
  const int num_blocks = (num_snps + (block_size-1)) / block_size;
  std::vector<int> idx1(block_elems, 0), idx2(block_elems, 0);
  for (int i = 0; i < block_size; i++) {
    for (int j = 0; j < block_size; j++) {
      idx1[i*block_size + j] = i;
      idx2[i*block_size + j] = j;
    }
  }

  PosixFile bedfile(bfile + ".bed", "rb");

  LdMatrixCsrChunk ld_matrix_csr_chunk;
  ld_matrix_csr_chunk.snp_index_from_inclusive_ = 0;
  ld_matrix_csr_chunk.snp_index_to_exclusive_ = num_snps;
  ld_matrix_csr_chunk.chr_label_ = 0;

  std::valarray<float> ld_tag_r2_sum(0.0, num_snps);
  std::valarray<float> ld_tag_r2_sum_adjust_for_hvec(0.0, num_snps);
  std::vector<float> freqvec(num_snps, 0.0);

  PlinkLdBedFileChunk chunk_fixed, chunk_var, *chunk_var_ptr;
  for (int block_idx = 0; block_idx < num_blocks; block_idx++) {
    const int block_istart = block_idx * block_size;
    const int block_iend = std::min(block_istart + block_size, num_snps);
    const int block_isize = block_iend - block_istart;
    if (0 != chunk_fixed.init(num_subj, block_istart, block_isize, bedfile.handle()))
      BGMG_THROW_EXCEPTION(::std::runtime_error("error while reading .bed file"));

    // save allele frequencies
    for (int block_snp_index = 0; block_snp_index < block_isize; block_snp_index++)
      freqvec[block_istart + block_snp_index] = chunk_fixed.freq()[block_snp_index];

    for (int block_jdx = block_idx; block_jdx < num_blocks; block_jdx++) {
    LOG << " processing block " << (block_idx+1) << "x" << (block_jdx+1) << " of " << num_blocks << "x" << num_blocks << "... ";

      const int block_jstart = block_jdx * block_size;
      const int block_jend = std::min(block_jstart + block_size, num_snps);
      const int block_jsize = block_jend - block_jstart;
      if (block_jdx == block_idx) {
        chunk_var_ptr = &chunk_fixed;  // reuse the block
      } else {
        if (0 != chunk_var.init(num_subj, block_jstart, block_jsize, bedfile.handle()))
          BGMG_THROW_EXCEPTION(::std::runtime_error("error while reading .bed file"));
        chunk_var_ptr = &chunk_var;
      }

      // There are many alternatives of collecting tail LD scores:
      // - collecting raw r2, versus r2 adjusted for heterozygosity
      // - collecting withing a range r2min < r2 < r2max;
      // - collecting raw r2, versus unbiased r2 estimates 
      // - collect r2 only when they pass significance at current sample size
      // For now, we choose to collect both raw r2 and r2 adjusted for heterozygosity,
      // within a range ldscore_r2min <= r2 < r2max,
      // without correcting for the bias or filtering based on p-value.
#pragma omp parallel
      {
        std::vector<std::tuple<int, int, packed_r_value>> local_coo_ld; // snp, tag, r2
        std::valarray<float> local_ld_tag_r2_sum(0.0, num_snps);
        std::valarray<float> local_ld_tag_r2_sum_adjust_for_hvec(0.0, num_snps);

#pragma omp for schedule(static)
        for (int k = 0; k < block_elems; k++) {
          int block_snp_index = k / block_size;
          int block_snp_jndex = k % block_size;
          int global_snp_index = block_snp_index + block_istart;
          int global_snp_jndex = block_snp_jndex + block_jstart;
          if (global_snp_jndex <= global_snp_index || global_snp_jndex >= num_snps) continue;
          float ld_corr = (float)PlinkLdBedFileChunk::calculate_ld_corr(chunk_fixed, *chunk_var_ptr, block_snp_index, block_snp_jndex);
          float ld_r2 = ld_corr * ld_corr;
          float ld_r4 = ld_r2 * ld_r2;

          if ((ldscore_r2min <= ld_r2) && (ld_r2 < r2_min)) {
            const float hval_at_index = chunk_fixed.hetval(block_snp_index);
            const float hval_at_jndex = chunk_var_ptr->hetval(block_snp_jndex);
            local_ld_tag_r2_sum[global_snp_index] += ld_r2;  // note that i-th SNP is adjusted for het of j-th SNP
            local_ld_tag_r2_sum[global_snp_jndex] += ld_r2;  // and vice versa.
            local_ld_tag_r2_sum_adjust_for_hvec[global_snp_index] += ld_r2 * hval_at_jndex;
            local_ld_tag_r2_sum_adjust_for_hvec[global_snp_jndex] += ld_r2 * hval_at_index;
          }
          if (ld_r2 >= r2_min) {
            local_coo_ld.push_back(std::make_tuple(global_snp_index, global_snp_jndex, ld_corr));
          }
        }
#pragma omp critical 
        {
          ld_tag_r2_sum += local_ld_tag_r2_sum;
          ld_tag_r2_sum_adjust_for_hvec += local_ld_tag_r2_sum_adjust_for_hvec;
          ld_matrix_csr_chunk.coo_ld_.insert( ld_matrix_csr_chunk.coo_ld_.end(), local_coo_ld.begin(), local_coo_ld.end() );
        }
      }
    }
  }

  ld_matrix_csr_chunk.set_ld_r2_csr(nullptr); 
  std::vector<float> ld_tag_r2_sum_adjust_for_hvec_vec, ld_tag_r2_sum_vec;
  ld_tag_r2_sum_vec.assign(std::begin(ld_tag_r2_sum), std::end(ld_tag_r2_sum));
  ld_tag_r2_sum_adjust_for_hvec_vec.assign(std::begin(ld_tag_r2_sum_adjust_for_hvec), std::end(ld_tag_r2_sum_adjust_for_hvec));

  save_ld_matrix(ld_matrix_csr_chunk, freqvec, ld_tag_r2_sum_vec, ld_tag_r2_sum_adjust_for_hvec_vec, outfile);

  LOG << ">generate_ld_matrix_from_bed_file(bfile=" << bfile << ", r2_min=" << r2_min << ", ldscore_r2min=" << ldscore_r2min << "), nnz=" << ld_matrix_csr_chunk.csr_ld_r_.size() <<", elapsed time " << timer.elapsed_ms() << "ms";
}

// reader must know the type
template<typename T>
void save_vector(std::ofstream& os, const std::vector<T>& vec) {
  size_t numel = vec.size();
  os.write(reinterpret_cast<const char*>(&numel), sizeof(size_t));
  os.write(reinterpret_cast<const char*>(&vec[0]), numel * sizeof(T));
}

template<typename T>
void save_value(std::ofstream& os, T value) {
  os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
void load_vector(std::ifstream& is, std::vector<T>* vec) {
  size_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(size_t));
  vec->resize(numel);
  is.read(reinterpret_cast<char*>(&(*vec)[0]), numel * sizeof(T));
}

template<typename T>
void load_value(std::ifstream& is, T* value) {
  is.read(reinterpret_cast<char*>(value), sizeof(T));
}

void save_ld_matrix(const LdMatrixCsrChunk& chunk,
                    const std::vector<float>& freqvec,
                    const std::vector<float>& ld_tag_r2_sum,
                    const std::vector<float>& ld_tag_r2_sum_adjust_for_hvec,
                    std::string filename) {
  std::ofstream os(filename, std::ofstream::binary);
  if (!os) BGMG_THROW_EXCEPTION(std::runtime_error(::std::runtime_error("can't open" + filename)));

  LOG << ">save_ld_matrix(filename=" << filename << "), format version " << LD_MATRIX_FORMAT_VERSION;

  size_t format_version = LD_MATRIX_FORMAT_VERSION;
  os.write(reinterpret_cast<const char*>(&format_version), sizeof(format_version));

  save_value(os, chunk.snp_index_from_inclusive_);
  save_value(os, chunk.snp_index_to_exclusive_);
  save_vector(os, chunk.csr_ld_snp_index_);
  save_vector(os, chunk.csr_ld_tag_index_offset_);
  save_vector(os, chunk.csr_ld_tag_index_packed_);
  save_vector(os, chunk.csr_ld_r_);

  save_vector(os, freqvec);
  save_vector(os, ld_tag_r2_sum);
  save_vector(os, ld_tag_r2_sum_adjust_for_hvec);

  os.close();

  LOG << "<save_ld_matrix(filename=" << filename << ")...";
}

// load LD matrix from an old format
void load_ld_matrix_version0(std::string filename,
                             std::vector<int>* snp_index,
                             std::vector<int>* tag_index,
                             std::vector<float>* r2) {
  LOG << ">load_ld_matrix_version0(filename=" << filename << ")";

  std::ifstream is(filename, std::ifstream::binary);
  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't open" + filename));
  if (sizeof(int) != 4) BGMG_THROW_EXCEPTION("sizeof(int) != 4, internal error in BGMG cpp"); // int -> int32_t

  int64_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(int64_t));
  LOG << " load_ld_matrix_version0(filename=" << filename << "), reading " << numel << " elements...";

  snp_index->resize(numel);
  tag_index->resize(numel);
  r2->resize(numel);

  is.read(reinterpret_cast<char*>(&snp_index->at(0)), numel * sizeof(int));
  is.read(reinterpret_cast<char*>(&tag_index->at(0)), numel * sizeof(int));
  is.read(reinterpret_cast<char*>(&r2->at(0)), numel * sizeof(float));

  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't read from " + filename));
  is.close();

  LOG << "<load_ld_matrix_version0(filename=" << filename << ")";
}

void load_ld_matrix(std::string filename,
                    LdMatrixCsrChunk* chunk,
                    std::vector<float>* freqvec,
                    std::vector<float>* ld_tag_r2_sum,
                    std::vector<float>* ld_tag_r2_sum_adjust_for_hvec) {
  LOG << ">load_ld_matrix(filename=" << filename << ")";

  std::ifstream is(filename, std::ifstream::binary);
  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't open" + filename));

  size_t format_version;
  is.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
  if (format_version <= 0 || format_version > LD_MATRIX_FORMAT_VERSION) throw("Unable to read an old format version");

  load_value(is, &chunk->snp_index_from_inclusive_);
  load_value(is, &chunk->snp_index_to_exclusive_);
  load_vector(is, &chunk->csr_ld_snp_index_);
  load_vector(is, &chunk->csr_ld_tag_index_offset_);
  load_vector(is, &chunk->csr_ld_tag_index_packed_);
  load_vector(is, &chunk->csr_ld_r_);

  load_vector(is, freqvec);
  load_vector(is, ld_tag_r2_sum);
  load_vector(is, ld_tag_r2_sum_adjust_for_hvec);

  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't read from " + filename));
  is.close();

  LOG << "<load_ld_matrix(filename=" << filename << "), format version " << format_version;
}
