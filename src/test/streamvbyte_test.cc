#include "gtest\gtest.h"

#include "streamvbyte.h"
#include "streamvbytedelta.h"

#include "TurboPFor/vsimple.h"

#include <vector>

namespace {

// size_t streamvbyte_encode(const uint32_t *in, uint32_t length, uint8_t *out);
// static size_t streamvbyte_max_compressedbytes(const uint32_t length) {
// size_t streamvbyte_decode(const uint8_t *in, uint32_t *out, uint32_t length);
// size_t streamvbyte_delta_encode(const uint32_t *in, uint32_t length, uint8_t *out, uint32_t prev);
// size_t streamvbyte_delta_decode(const uint8_t *in, uint32_t *out, uint32_t length, uint32_t prev);
//  unsigned char *vsenc32(unsigned       *__restrict in, size_t n, unsigned char  *__restrict out);
//  unsigned char *vsdec32(unsigned char  *__restrict in, size_t n, unsigned       *__restrict out);

void TurboPForVSimpleTest(int max_size) {
  for (int seq = 0; seq < 5; seq++) {
    std::vector<uint32_t> data;
    int numel = (16 + (rand() % 512));
    for (int i = 0; i < numel; i++)
      data.push_back(rand() % max_size);

    size_t buflen = 10000 + data.size() * sizeof(uint32_t);
    std::vector<unsigned char> buffer(buflen, 0);

    uint32_t prev = 0;
    unsigned char *outptr = vsenc32(&data[0], data.size(), &buffer[0]);
    size_t encoded_bytes = outptr - &buffer[0];
    ASSERT_LE(encoded_bytes, buflen);

    std::vector<uint32_t> data2(data.size(), 0);
    outptr = vsdec32(&buffer[0], data.size(), &data2[0]);

    size_t decoded_bytes = outptr - &buffer[0];
    ASSERT_EQ(decoded_bytes, encoded_bytes);

    for (int i = 0; i < data.size(); i++) ASSERT_EQ(data[i], data2[i]);
    printf("len %i, max values %i, ratio = %.3f\n", data.size(), max_size, float(encoded_bytes) / (float)(data.size() * sizeof(uint32_t)));
  }
}

void StreamVByteTest(bool use_delta, int max_size) {
  std::vector<uint32_t> data;
  uint32_t val = rand();
  for (int seq = 0; seq < 10; seq++)
    for (int i = 0; i < 1000; i++, val += rand() % max_size)
      data.push_back(val);

  size_t buflen = streamvbyte_max_compressedbytes(data.size());
  std::vector<uint8_t> buffer(buflen, 0);

  uint32_t prev = 0;
  size_t encoded_bytes = (use_delta ? 
    streamvbyte_delta_encode(&data[0], data.size(), &buffer[0], prev) : 
    streamvbyte_encode(&data[0], data.size(), &buffer[0]));

  std::vector<uint32_t> data2(data.size(), 0);
  size_t decoded_bytes = (use_delta ?
    streamvbyte_delta_decode(&buffer[0], &data2[0], data.size(), prev) : 
    streamvbyte_decode(&buffer[0], &data2[0], data.size()));

  ASSERT_EQ(decoded_bytes, encoded_bytes);
  for (int i = 0; i < data.size(); i++) ASSERT_EQ(data[i], data2[i]);
  printf("max values %i, ratio = %.3f\n", max_size, float(encoded_bytes) / ((float)data.size() * sizeof(uint32_t)));
}


// bgmg-test.exe --gtest_filter=Compression.StreamVByteEncodeTest
TEST(Compression, StreamVByteEncodeTest) {
  for (int i=4; i<=2048; i*=2)
    StreamVByteTest(false, i);
}

// bgmg-test.exe --gtest_filter=Compression.StreamVByteDeltaEncodeTest
TEST(Compression, StreamVByteDeltaEncodeTest) {
  for (int i = 4; i <= 2048; i *= 2)
    StreamVByteTest(true, i);
}

// bgmg-test.exe --gtest_filter=Compression.TurboPForVSimpleTest
TEST(Compression, TurboPForVSimpleTest) {
  for (int i = 4; i <= 2048; i *= 2)
    TurboPForVSimpleTest(i);
}

}  // namespace