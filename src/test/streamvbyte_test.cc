#include "gtest/gtest.h"

#include "TurboPFor/vsimple.h"

#include <vector>

namespace {

//  unsigned char *vsenc32(unsigned       *__restrict in, size_t n, unsigned char  *__restrict out);
//  unsigned char *vsdec32(unsigned char  *__restrict in, size_t n, unsigned       *__restrict out);

#define ROUND_UP(_n_, _a_) (((_n_) + ((_a_)-1)) & ~((_a_)-1))
#define P4NENC_BOUND(n, size) ((n + 127) / 128 + (n + 32) * (size))
#define P4NDEC_BOUND(n, size) (ROUND_UP(n, 32) * (size))

void TurboPForVSimpleTest(int max_size) {
  printf("Max val: %i, ", max_size);
  for (int seq = 0; seq < 5; seq++) {
    std::vector<uint32_t> data;
    int numel = (16 + (rand() % 512));
    for (int i = 0; i < numel; i++)
      data.push_back(rand() % max_size);

    size_t buflen = P4NENC_BOUND(data.size(), sizeof(uint32_t));
    std::vector<unsigned char> buffer(buflen, 0);

    uint32_t prev = 0;
    unsigned char *outptr = vsenc32(&data[0], data.size(), &buffer[0]);
    size_t encoded_bytes = outptr - &buffer[0];
    ASSERT_LE(encoded_bytes, buflen);

    std::vector<uint32_t> data2(32 + ROUND_UP(data.size(), 32), 0);
    outptr = vsdec32(&buffer[0], data.size(), &data2[0]);

    size_t decoded_bytes = outptr - &buffer[0];
    ASSERT_EQ(decoded_bytes, encoded_bytes);

    for (int i = 0; i < data.size(); i++) ASSERT_EQ(data[i], data2[i]);
    printf("%i@%.1f%; ", data.size(), 100 * float(encoded_bytes) / (float)(data.size() * sizeof(uint32_t)));
  }
  printf("\n");
}

// bgmg-test.exe --gtest_filter=Compression.TurboPForVSimpleTest
TEST(Compression, TurboPForVSimpleTest) {
  for (int i = 4; i <= 2048; i *= 2)
    TurboPForVSimpleTest(i);
}

}  // namespace
