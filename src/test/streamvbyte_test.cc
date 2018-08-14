#include "gtest\gtest.h"

#include "streamvbyte.h"
#include "streamvbytedelta.h"

#include <vector>

namespace {

// size_t streamvbyte_encode(const uint32_t *in, uint32_t length, uint8_t *out);
// static size_t streamvbyte_max_compressedbytes(const uint32_t length) {
// size_t streamvbyte_decode(const uint8_t *in, uint32_t *out, uint32_t length);
// size_t streamvbyte_delta_encode(const uint32_t *in, uint32_t length, uint8_t *out, uint32_t prev);
// size_t streamvbyte_delta_decode(const uint8_t *in, uint32_t *out, uint32_t length, uint32_t prev);

void StreamVByteTest(bool use_delta) {
  std::vector<uint32_t> data;
  uint32_t val = rand();
  for (int seq = 0; seq < 10; seq++)
    for (int i = 0; i < 1000; i++, val += rand() % 512)
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
}

// bgmg-test.exe --gtest_filter=StreamVByte.EncodeTest
TEST(StreamVByte, EncodeTest) {
  StreamVByteTest(false);
}

// bgmg-test.exe --gtest_filter=StreamVByte.DeltaEncodeTest
TEST(StreamVByte, DeltaEncodeTest) {
  StreamVByteTest(true);
}

}  // namespace