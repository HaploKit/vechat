// Copyright (c) 2020 Robert Vaser

#include "bioparser/parser.hpp"

#include "biosoup/sequence.hpp"
#include "gtest/gtest.h"

#include "bioparser/fasta_parser.hpp"

namespace bioparser {
namespace test {

TEST(BioparserParserTest, Create) {
  try {
    auto p = Parser<biosoup::Sequence>::Create<FastaParser>("");
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::Parser::Create] error: unable to open file ");
  }
}

}  // namespace test
}  // namespace bioparser
