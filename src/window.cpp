/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>

#include "window.hpp"

#include "spoa/spoa.hpp"

namespace racon
{

    std::shared_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
                                         const char *backbone, uint32_t backbone_length, const char *quality,
                                         uint32_t quality_length)
    {

        if (backbone_length == 0 || backbone_length != quality_length)
        {
            fprintf(stderr, "[racon::createWindow] error: "
                            "empty backbone sequence/unequal quality length!\n");
            exit(1);
        }

        return std::shared_ptr<Window>(new Window(id, rank, type, backbone,
                                                  backbone_length, quality, quality_length));
    }

    Window::Window(uint64_t id, uint32_t rank, WindowType type, const char *backbone,
                   uint32_t backbone_length, const char *quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), sequences_(),
          qualities_(), positions_()
    {

        sequences_.emplace_back(backbone, backbone_length);
        qualities_.emplace_back(quality, quality_length);
        positions_.emplace_back(0, 0);
    }

    Window::~Window()
    {
    }

    void Window::add_layer(const char *sequence, uint32_t sequence_length,
                           const char *quality, uint32_t quality_length, uint32_t begin, uint32_t end)
    {

        if (sequence_length == 0 || begin == end)
        {
            return;
        }

        if (quality != nullptr && sequence_length != quality_length)
        {
            fprintf(stderr, "[racon::Window::add_layer] error: "
                            "unequal quality size!\n");
            exit(1);
        }
        if (begin >= end || begin > sequences_.front().second || end > sequences_.front().second)
        {
            fprintf(stderr, "[racon::Window::add_layer] error: "
                            "layer begin and end positions are invalid!\n");
            exit(1);
        }

        sequences_.emplace_back(sequence, sequence_length);
        qualities_.emplace_back(quality, quality_length);
        positions_.emplace_back(begin, end);
    }

    bool Window::generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine,
                                    bool trim)
    {

        if (sequences_.size() < 3)
        {
            consensus_ = std::string(sequences_.front().first, sequences_.front().second);
            return false;
        }

        spoa::Graph graph{};
        graph.AddAlignment(
            spoa::Alignment(),
            sequences_.front().first, sequences_.front().second,
            qualities_.front().first, qualities_.front().second);

        std::vector<uint32_t> rank;
        rank.reserve(sequences_.size());
        for (uint32_t i = 0; i < sequences_.size(); ++i)
        {
            rank.emplace_back(i);
        }

        std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) { return positions_[lhs].first < positions_[rhs].first; });

        std::cerr << "Nohap sequences_[0]:" << sequences_.front().first << std::endl;
        std::cerr << "Nohap sequences_[rank[0]]:" << sequences_[rank[0]].first << std::endl;

        uint32_t offset = 0.01 * sequences_.front().second;
        for (uint32_t j = 1; j < sequences_.size(); ++j) //j starts from 1, the 0th is the backbone
        {
            uint32_t i = rank[j];

            spoa::Alignment alignment;
            if (positions_[i].first < offset && positions_[i].second >
                                                    sequences_.front().second - offset)
            {
                alignment = alignment_engine->Align(
                    sequences_[i].first, sequences_[i].second, graph);
            }
            else
            {
                // I guess this is used when the sequence to be aligned is short (NGS or ends of long read)
                // such that it is not good to perform global alignment on the whole POA graph,
                // but it is fine if only on the subgraph covered by the target sequence.
                std::vector<const spoa::Graph::Node *> mapping;
                auto subgraph = graph.Subgraph(
                    positions_[i].first,
                    positions_[i].second,
                    &mapping);
                alignment = alignment_engine->Align(
                    sequences_[i].first, sequences_[i].second, subgraph);
                subgraph.UpdateAlignment(mapping, &alignment);
            }

            if (qualities_[i].first == nullptr)
            {
                graph.AddAlignment(
                    alignment,
                    sequences_[i].first, sequences_[i].second);
            }
            else
            {
                graph.AddAlignment(
                    alignment,
                    sequences_[i].first, sequences_[i].second,
                    qualities_[i].first, qualities_[i].second);
            }
        }

        std::vector<uint32_t> coverages;
        consensus_ = graph.GenerateConsensus(&coverages);
        std::cerr << "consensus_ len:" << consensus_.length() << std::endl;

        if (type_ == WindowType::kTGS && trim)
        {
            uint32_t average_coverage = (sequences_.size() - 1) / 2; //maybe too strict in our case

            int32_t begin = 0, end = consensus_.size() - 1;
            for (; begin < static_cast<int32_t>(consensus_.size()); ++begin)
            {
                if (coverages[begin] >= average_coverage)
                {
                    break;
                }
            }
            for (; end >= 0; --end)
            {
                if (coverages[end] >= average_coverage)
                {
                    break;
                }
            }

            if (begin >= end)
            {
                fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                                "contig %lu might be chimeric in window %u!\n",
                        id_, rank_);
            }
            else
            {
                consensus_ = consensus_.substr(begin, end - begin + 1);
            }
        }

        return true;
    }

    bool Window::generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine,
                                    bool trim, bool haplotype)
    {
        std::cerr << "Using --haplotype for error correction !!! " << std::endl;
        //function overloading
        if (!haplotype)
        {
            throw std::invalid_argument(
                "[spoa::Window::generate_consensus] error: "
                "invalid haplotype mode!");
        }
        if (sequences_.size() < 3)
        {
            consensus_ = std::string(sequences_.front().first, sequences_.front().second);
            return false;
        }

        // it seems that sequences_ stores the whole sequence of read, not window sequence
        spoa::Graph graph{};
        graph.AddAlignment(
            spoa::Alignment(),
            sequences_.front().first, sequences_.front().second,
            qualities_.front().first, qualities_.front().second);

        std::vector<uint32_t> rank;
        rank.reserve(sequences_.size());
        for (uint32_t i = 0; i < sequences_.size(); ++i)
        {
            rank.emplace_back(i);
        }

        std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) { return positions_[lhs].first < positions_[rhs].first; });

        uint32_t offset = 0.01 * sequences_.front().second;
        // the original POA graph construction
        for (uint32_t j = 1; j < sequences_.size(); ++j)
        {
            uint32_t i = rank[j];

            spoa::Alignment alignment;
            if (positions_[i].first < offset && positions_[i].second >
                                                    sequences_.front().second - offset)
            {
                alignment = alignment_engine->Align(
                    sequences_[i].first, sequences_[i].second,
                    graph);
            }
            else
            {
                // I guess this is used when the sequence to be aligned is short (NGS or ends of long read)
                // such that it is not good to perform global alignment on the whole POA graph,
                // but it is fine if only on the subgraph covered by the target sequence.
                std::vector<const spoa::Graph::Node *> mapping;
                auto subgraph = graph.Subgraph(
                    positions_[i].first,
                    positions_[i].second,
                    &mapping);
                alignment = alignment_engine->Align(
                    sequences_[i].first, sequences_[i].second,
                    subgraph);
                //replace the node IDs(first column, subgraph) in the alignment
                // with the node IDs in the original graph
                subgraph.UpdateAlignment(mapping, &alignment);
            }

            if (qualities_[i].first == nullptr)
            {
                graph.AddAlignment(
                    alignment,
                    sequences_[i].first, sequences_[i].second);
            }
            else
            {
                graph.AddAlignment(
                    alignment,
                    sequences_[i].first, sequences_[i].second,
                    qualities_[i].first, qualities_[i].second);
            }
        }

        std::cerr << "old sequences_[0]:" << sequences_[0].first << std::endl;
        std::cerr << "old sequences_[rank[0]]:" << sequences_[rank[0]].first << std::endl;

        // start to prune the graph

        int64_t min_weight = 3;
        double min_confidence = 0.05;
        double min_support = 0.1;
        std::uint32_t num_prune = 2;

        std::cerr << "Pruning graph " << 1 << "th...\n";
        graph.PruneGraph(min_weight, min_confidence, min_support);

        spoa::Graph lagestsubgraph{};
        lagestsubgraph = graph.LargestSubgraph();
        graph.Clear();

        spoa::Graph *p;
        p = &lagestsubgraph;

        //TODO, try to consider SubGraph align !!

        // for (std::uint32_t k = 0; k < num_prune - 1; ++k)
        // {
            // std::cerr << "Pruning graph " << k + 2 << "th...\n";

            for (uint32_t j = 1; j < sequences_.size(); ++j)
            {
                uint32_t i = rank[j];

                spoa::Alignment alignment;
                alignment = alignment_engine->Align(sequences_[i].first, sequences_[i].second, *p);
                // std::cerr << "alignment size1: " << alignment.size()<<std::endl;

                std::vector<std::uint32_t> weights;
                for (std::uint32_t n = 0; n < sequences_[i].second; ++n)
                {
                    weights.emplace_back(1); // TODO, consider quality socre
                }
                (*p).AddWeights(alignment, sequences_[i].first, sequences_[i].second, weights);
            }

            (*p).PruneGraph(min_weight, min_confidence, min_support);

            spoa::Graph lagestsubgraph2{};
            lagestsubgraph2 = (*p).LargestSubgraph();
            // lagestsubgraph.Clear();
            p = &lagestsubgraph2;
            std::cerr << "lagestsubgraph2 node size:" << (*p).nodes().size() << std::endl;
        // }

        std::cerr << "lagestsubgraph node size:" << (*p).nodes().size() << std::endl;
        std::cerr << "lagestsubgraph edge size:" << (*p).edges().size() << std::endl;

        //only need to correct the target read (the first one, I guess)
        std::cerr << "sequences_[0]:" << sequences_[0].first << std::endl;
        std::cerr << "sequences_[rank[0]]:" << sequences_[rank[0]].first << std::endl;

        auto alignment = alignment_engine->Align(sequences_.front().first, sequences_.front().second, *p);
        // auto alignment = alignment_engine->Align(sequences_[rank[0]].first, sequences_[rank[0]].second, *p);

        consensus_ = (*p).GenerateCorrectedSequence(alignment);
        std::cerr << "alignment size2: " << alignment.size() << std::endl;
        // for (const auto &it_align : alignment)
        // {
        //     std::cerr << "alignment: " << it_align.first << "\t" << it_align.second << std::endl;
        // }
        std::cerr << ">Window consensus: " << consensus_ << std::endl;

        /*  TODO: trim consensus based on base coverages to avoid chimeric sequences    

    std::vector<uint32_t> coverages;
    consensus_ = graph.GenerateConsensus(&coverages);

    if (type_ == WindowType::kTGS && trim) {
        uint32_t average_coverage = (sequences_.size() - 1) / 2; //maybe too strict in our case

        int32_t begin = 0, end = consensus_.size() - 1;
        for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
            if (coverages[begin] >= average_coverage) {
                break;
            }
        }
        for (; end >= 0; --end) {
            if (coverages[end] >= average_coverage) {
                break;
            }
        }

        if (begin >= end) {
            fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                "contig %lu might be chimeric in window %u!\n", id_, rank_);
        } else {
            consensus_ = consensus_.substr(begin, end - begin + 1);
        }
    } 
*/

        return true;
    }
}
