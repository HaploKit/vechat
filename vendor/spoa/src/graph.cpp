// Copyright (c) 2020 Robert Vaser

#include "spoa/graph.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <stack>
#include <stdexcept>
#include <unordered_set>

namespace spoa
{

  Graph::Node::Node(std::uint32_t id, std::uint32_t code)
      : id(id),
        code(code),
        inedges(),
        outedges(),
        aligned_nodes()
  {
  }

  Graph::Node *Graph::Node::Successor(std::uint32_t label) const
  {
    for (const auto &it : outedges)
    {
      auto jt = std::find(it->labels.begin(), it->labels.end(), label);
      if (jt != it->labels.end())
      {
        return it->head;
      }
    }
    return nullptr;
  }

  std::uint32_t Graph::Node::Coverage() const
  {
    std::unordered_set<std::uint32_t> labels;
    for (const auto &it : inedges)
    {
      std::copy(
          it->labels.begin(),
          it->labels.end(),
          std::inserter(labels, labels.end()));
    }
    for (const auto &it : outedges)
    {
      std::copy(
          it->labels.begin(),
          it->labels.end(),
          std::inserter(labels, labels.end()));
    }
    return labels.size();
  }

  Graph::Edge::Edge(
      Node *tail,
      Node *head,
      std::uint32_t label,
      std::uint32_t weight)
      : tail(tail),
        head(head),
        labels(1, label),
        weight(weight)
  {
  }

  void Graph::Edge::AddSequence(std::uint32_t label, std::uint32_t w)
  {
    labels.emplace_back(label);
    weight += w;
  }

  Graph::Graph()
      : num_codes_(0),
        coder_(256, -1),
        decoder_(256, -1),
        sequences_(),
        nodes_(),
        edges_(),
        rank_to_node_(),
        consensus_()
  {
  }

  Graph::Node *Graph::AddNode(std::uint32_t code)
  {
    nodes_.emplace_back(new Node(nodes_.size(), code));
    return nodes_.back().get();
  }

  void Graph::AddEdge(Node *tail, Node *head, std::uint32_t weight)
  {
    for (const auto &it : tail->outedges)
    {
      if (it->head == head)
      {
        it->AddSequence(sequences_.size(), weight);
        return;
      }
    }
    edges_.emplace_back(new Edge(tail, head, sequences_.size(), weight));
    tail->outedges.emplace_back(edges_.back().get());
    head->inedges.emplace_back(edges_.back().get());
  }

  Graph::Node *Graph::AddSequence(
      const char *sequence,
      const std::vector<std::uint32_t> &weights,
      std::uint32_t begin,
      std::uint32_t end)
  {
    if (begin == end)
    {
      return nullptr;
    }
    Node *prev = nullptr;
    for (std::uint32_t i = begin; i < end; ++i)
    {
      auto curr = AddNode(coder_[sequence[i]]);
      if (prev)
      { // both nodes contribute to the weight
        AddEdge(prev, curr, weights[i - 1] + weights[i]);
      }
      prev = curr;
    }
    return nodes_[nodes_.size() - (end - begin)].get();
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const std::string &sequence,
      std::uint32_t weight)
  {
    AddAlignment(alignment, sequence.c_str(), sequence.size(), weight);
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const char *sequence, std::uint32_t sequence_len,
      std::uint32_t weight)
  {
    std::vector<std::uint32_t> weights(sequence_len, weight);
    AddAlignment(alignment, sequence, sequence_len, weights);
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const std::string &sequence,
      const std::string &quality)
  {
    AddAlignment(
        alignment,
        sequence.c_str(), sequence.size(),
        quality.c_str(), quality.size());
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const char *sequence, std::uint32_t sequence_len,
      const char *quality, std::uint32_t quality_len)
  {
    std::vector<std::uint32_t> weights;
    for (std::uint32_t i = 0; i < quality_len; ++i)
    {
      weights.emplace_back(quality[i] - 33); // Phred quality
    }
    AddAlignment(alignment, sequence, sequence_len, weights);
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const std::string &sequence,
      const std::vector<std::uint32_t> &weights)
  {
    AddAlignment(alignment, sequence.c_str(), sequence.size(), weights);
  }

  void Graph::AddAlignment(
      const Alignment &alignment,
      const char *sequence, std::uint32_t sequence_len,
      const std::vector<std::uint32_t> &weights)
  {
    if (sequence_len == 0)
    {
      return;
    }
    if (sequence_len != weights.size())
    {
      throw std::invalid_argument(
          "[spoa::Graph::AddAlignment] error: "
          "sequence and weights are of unequal size!");
    }

    for (std::uint32_t i = 0; i < sequence_len; ++i)
    {
      if (coder_[sequence[i]] == -1)
      {
        coder_[sequence[i]] = num_codes_;
        decoder_[num_codes_++] = sequence[i];
        std::cerr << "sequence[i]:" << sequence[i] << "\ti=" << i << "\tseq size:" << std::strlen(sequence) << std::endl;
      }
    }

    if (alignment.empty())
    {
      sequences_.emplace_back(AddSequence(sequence, weights, 0, sequence_len));
      TopologicalSort();
      return;
    }

    std::vector<std::uint32_t> valid;
    for (const auto &it : alignment)
    {
      if (it.second != -1)
      {
        if (it.second < 0 || it.second >= static_cast<std::int32_t>(sequence_len))
        { // NOLINT
          throw std::invalid_argument(
              "[spoa::Graph::AddAlignment] error: invalid alignment");
        }
        valid.emplace_back(it.second);
      }
    }
    if (valid.empty())
    {
      throw std::invalid_argument(
          "[spoa::Graph::AddAlignment] error: missing sequence in alignment");
    }

    // add unaligned bases
    Node *begin = AddSequence(sequence, weights, 0, valid.front());
    Node *prev = begin ? nodes_.back().get() : nullptr;
    Node *last = AddSequence(sequence, weights, valid.back() + 1, sequence_len);

    // add aligned bases
    for (const auto &it : alignment)
    {
      if (it.second == -1)
      {
        continue;
      }

      std::uint32_t code = coder_[sequence[it.second]];
      Node *curr = nullptr;
      if (it.first == -1)
      {
        curr = AddNode(code);
      }
      else
      {
        auto jt = nodes_[it.first].get();
        if (jt->code == code)
        {
          curr = jt;
        }
        else
        {
          for (const auto &kt : jt->aligned_nodes)
          {
            if (kt->code == code)
            {
              curr = kt;
              break;
            }
          }
          if (!curr)
          {
            curr = AddNode(code);
            for (const auto &kt : jt->aligned_nodes)
            {
              kt->aligned_nodes.emplace_back(curr);
              curr->aligned_nodes.emplace_back(kt);
            }
            jt->aligned_nodes.emplace_back(curr);
            curr->aligned_nodes.emplace_back(jt);
          }
        }
      }
      if (!begin)
      {
        begin = curr;
      }
      if (prev)
      { // both nodes contribute to weight
        AddEdge(prev, curr, weights[it.second - 1] + weights[it.second]);
      }
      prev = curr;
    }
    if (last)
    {
      AddEdge(prev, last, weights[valid.back()] + weights[valid.back() + 1]);
    }
    sequences_.emplace_back(begin);

    TopologicalSort();
  }

  void Graph::TopologicalSort()
  {
    rank_to_node_.clear();

    std::vector<std::uint8_t> marks(nodes_.size(), 0);
    std::vector<bool> ignored(nodes_.size(), 0);

    std::stack<Node *> stack;
    for (const auto &it : nodes_)
    {
      if (marks[it->id] != 0)
      {
        continue;
      }
      stack.push(it.get());
      while (!stack.empty())
      {
        auto curr = stack.top();
        bool is_valid = true;
        if (marks[curr->id] != 2)
        {
          for (const auto &jt : curr->inedges)
          {
            if (marks[jt->tail->id] != 2)
            {
              stack.push(jt->tail);
              is_valid = false;
            }
          }
          if (!ignored[curr->id])
          {
            for (const auto &jt : curr->aligned_nodes)
            {
              if (marks[jt->id] != 2)
              {
                stack.push(jt);
                ignored[jt->id] = true;
                is_valid = false;
              }
            }
          }

          assert((is_valid || marks[curr->id] != 1) && "Graph is not a DAG");

          if (is_valid)
          {
            marks[curr->id] = 2;
            if (!ignored[curr->id])
            {
              rank_to_node_.emplace_back(curr);
              for (const auto &jt : curr->aligned_nodes)
              {
                rank_to_node_.emplace_back(jt);
              }
            }
          }
          else
          {
            marks[curr->id] = 1;
          }
        }

        if (is_valid)
        {
          stack.pop();
        }
      }
    }

    assert(IsTopologicallySorted() && "Graph is not topologically sorted");
  }

  bool Graph::IsTopologicallySorted() const
  {
    assert(nodes_.size() == rank_to_node_.size() && "Topological sort not called "); // NOLINT

    std::vector<bool> visited(nodes_.size(), 0);
    for (const auto &it : rank_to_node_)
    {
      for (const auto &jt : it->inedges)
      {
        if (!visited[jt->tail->id])
        {
          return false;
        }
      }
      visited[it->id] = 1;
    }

    return true;
  }

  std::vector<std::uint32_t> Graph::InitializeMultipleSequenceAlignment(
      std::uint32_t *row_size) const
  {
    std::vector<std::uint32_t> dst(nodes_.size());
    std::uint32_t j = 0;
    for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i, ++j)
    {
      auto it = rank_to_node_[i];
      dst[it->id] = j;
      for (const auto &jt : it->aligned_nodes)
      {
        dst[jt->id] = j;
        ++i;
      }
    }
    if (row_size)
    {
      *row_size = j;
    }
    return dst;
  }

  std::vector<std::string> Graph::GenerateMultipleSequenceAlignment(
      bool include_consensus)
  {
    std::uint32_t row_size = 0;
    auto node_id_to_column = InitializeMultipleSequenceAlignment(&row_size);

    std::vector<std::string> dst;
    for (std::uint32_t i = 0; i < sequences_.size(); ++i)
    {
      std::string row(row_size, '-');
      auto it = sequences_[i];
      while (true)
      {
        row[node_id_to_column[it->id]] = decoder_[it->code];
        if (!(it = it->Successor(i)))
        {
          break;
        }
      }
      dst.emplace_back(row);
    }
    if (include_consensus)
    {
      TraverseHeaviestBundle();
      std::string row(row_size, '-');
      for (const auto &it : consensus_)
      {
        row[node_id_to_column[it->id]] = decoder_[it->code];
      }
      dst.emplace_back(row);
    }

    return dst;
  }

  std::string Graph::GenerateConsensus()
  {
    TraverseHeaviestBundle();
    std::string dst{};
    for (const auto &it : consensus_)
    {
      dst += decoder_[it->code];
    }
    return dst;
  }

  std::string Graph::GenerateConsensus(
      std::vector<std::uint32_t> *summary,
      bool verbose)
  {
    if (!summary)
    {
      throw std::invalid_argument(
          "[spoa::Graph::GenerateConsensus] error: invalid ptr to summary");
    }

    auto dst = GenerateConsensus();

    summary->clear();
    if (!verbose)
    {
      for (const auto &it : consensus_)
      {
        summary->emplace_back(0);
        summary->back() += it->Coverage();
        for (const auto &jt : it->aligned_nodes)
        {
          summary->back() += jt->Coverage();
        }
      }
    }
    else
    {
      summary->resize((num_codes_ + 1) * consensus_.size(), 0);
      auto node_id_to_column = InitializeMultipleSequenceAlignment();

      for (std::uint32_t i = 0; i < sequences_.size(); ++i)
      {
        Node *it = sequences_[i];
        std::uint32_t c = 0, p, column = node_id_to_column[it->id];
        bool is_gap = false;
        while (true)
        {
          for (; c < consensus_.size(); ++c)
          {
            if (node_id_to_column[consensus_[c]->id] < column)
            {
              continue;
            }
            else
            {
              if (node_id_to_column[consensus_[c]->id] == column)
              {
                if (is_gap)
                {
                  for (std::uint32_t j = p + 1; j < c; ++j)
                  {
                    ++(*summary)[num_codes_ * consensus_.size() + j];
                  }
                }
                is_gap = true;
                p = c;
                ++(*summary)[it->code * consensus_.size() + c];
              }
              break;
            }
          }
          if (c == consensus_.size() || !(it = it->Successor(i)))
          {
            break;
          }
          column = node_id_to_column[it->id];
        }
      }
    }

    return dst;
  }

  void Graph::TraverseHeaviestBundle()
  {
    if (rank_to_node_.empty())
    {
      return;
    }

    std::vector<Node *> predecessors(nodes_.size(), nullptr);
    std::vector<std::int64_t> scores(nodes_.size(), -1);
    Node *max = nullptr;

    for (const auto &it : rank_to_node_)
    {
      for (const auto &jt : it->inedges)
      {
        if ((scores[it->id] < jt->weight) ||
            (scores[it->id] == jt->weight && scores[predecessors[it->id]->id] <= scores[jt->tail->id]))
        { // NOLINT
          scores[it->id] = jt->weight;
          predecessors[it->id] = jt->tail;
        }
      }
      if (predecessors[it->id])
      {
        scores[it->id] += scores[predecessors[it->id]->id];
      }
      if (!max || scores[max->id] < scores[it->id])
      {
        max = it;
      }
    }

    if (!max->outedges.empty())
    {
      std::vector<std::uint32_t> node_id_to_rank(nodes_.size(), 0);
      for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i)
      {
        node_id_to_rank[rank_to_node_[i]->id] = i;
      }
      while (!max->outedges.empty())
      {
        max = BranchCompletion(node_id_to_rank[max->id], &scores, &predecessors);
      }
    }

    // traceback
    consensus_.clear();
    while (predecessors[max->id])
    {
      consensus_.emplace_back(max);
      max = predecessors[max->id];
    }
    consensus_.emplace_back(max);
    std::reverse(consensus_.begin(), consensus_.end());
  }

  Graph::Node *Graph::BranchCompletion(
      std::uint32_t rank,
      std::vector<std::int64_t> *scores,
      std::vector<Node *> *predecessors)
  {
    auto start = rank_to_node_[rank];
    for (const auto &it : start->outedges)
    {
      for (const auto &jt : it->head->inedges)
      {
        if (jt->tail != start)
        {
          (*scores)[jt->tail->id] = -1;
        }
      }
    }

    Node *max = nullptr;
    for (std::uint32_t i = rank + 1; i < rank_to_node_.size(); ++i)
    {
      auto it = rank_to_node_[i];
      (*scores)[it->id] = -1;
      (*predecessors)[it->id] = nullptr;

      for (const auto &jt : it->inedges)
      {
        if ((*scores)[jt->tail->id] == -1)
        {
          continue;
        }
        if (((*scores)[it->id] < jt->weight) ||
            ((*scores)[it->id] == jt->weight && (*scores)[(*predecessors)[it->id]->id] <= (*scores)[jt->tail->id]))
        { // NOLINT
          (*scores)[it->id] = jt->weight;
          (*predecessors)[it->id] = jt->tail;
        }
      }
      if ((*predecessors)[it->id])
      {
        (*scores)[it->id] += (*scores)[(*predecessors)[it->id]->id];
      }
      if (!max || (*scores)[max->id] < (*scores)[it->id])
      {
        max = it;
      }
    }

    return max;
  }

  std::vector<bool> Graph::ExtractSubgraph(const Node *begin, const Node *end) const
  { // NOLINT
    std::vector<bool> dst(nodes_.size(), false);
    std::stack<const Node *> stack;
    stack.push(begin);
    int x = 0;
    while (!stack.empty())
    {
      auto curr = stack.top();

      stack.pop();
      /*       std::cerr << "x=" << x++ << std::endl;
      std::cerr << "curr->inedges.size " << curr->inedges.size() << std::endl; 
      std::cerr << "curr->code " << curr->code << std::endl; 
      std::cerr << "curr->id " << curr->id << std::endl; //bug, id=1413562689

      std::cerr << "end->id " << end->id << std::endl; */
      if (!dst[curr->id] && curr->id >= end->id)
      {
        // std::cerr << "fuck" << std::endl;
        for (const auto &it : curr->inedges)
        {
          // std::cerr << "fuck2" << std::endl;
          stack.push(it->tail);
        }
        for (const auto &it : curr->aligned_nodes)
        {
          stack.push(it);
        }
        dst[curr->id] = true;
      }
    }

    return dst;
  }

  Graph Graph::Subgraph(
      std::uint32_t begin,
      std::uint32_t end,
      std::vector<const Node *> *subgraph_to_graph) const
  {
    if (!subgraph_to_graph)
    {
      throw std::invalid_argument(
          "[spoa::Graph::Subgraph] error: invalid ptr to subgraph_to_graph");
    }

    // std::cerr << "subgraph func1..." << std::endl;
    // std::cerr << "end.." << nodes_[end].get() << std::endl;
    // std::cerr << "begin.." << nodes_[begin].get() << std::endl;
    // std::cerr << "subgraph func2.." << std::endl;
    auto is_in_subgraph = ExtractSubgraph(nodes_[end].get(), nodes_[begin].get());

    // init subgraph
    Graph subgraph{};
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;
    // subgraph.sequences_ = TODO(rvaser) maybe add sequences

    // create a map from subgraph nodes to graph nodes and vice versa
    subgraph_to_graph->clear();
    subgraph_to_graph->resize(nodes_.size(), nullptr);

    std::vector<Node *> graph_to_subgraph(nodes_.size(), nullptr);

    for (const auto &it : nodes_)
    {
      if (!is_in_subgraph[it->id])
      {
        continue;
      }
      subgraph.AddNode(it->code);
      graph_to_subgraph[it->id] = subgraph.nodes_.back().get();
      (*subgraph_to_graph)[subgraph.nodes_.back()->id] = it.get();
    }

    // connect nodes
    for (const auto &it : nodes_)
    {
      if (!is_in_subgraph[it->id])
      {
        continue;
      }
      auto jt = graph_to_subgraph[it->id];
      for (const auto &kt : it->inedges)
      {
        if (graph_to_subgraph[kt->tail->id])
        {
          subgraph.AddEdge(graph_to_subgraph[kt->tail->id], jt, kt->weight);
        }
      }
      for (const auto &kt : it->aligned_nodes)
      {
        if (graph_to_subgraph[kt->id])
        {
          jt->aligned_nodes.emplace_back(graph_to_subgraph[kt->id]);
        }
      }
    }

    subgraph.TopologicalSort();

    return subgraph;
  }

  void Graph::UpdateAlignment(
      const std::vector<const Node *> &subgraph_to_graph,
      Alignment *alignment) const
  {
    for (auto &it : *alignment)
    {
      if (it.first != -1)
      {
        it.first = subgraph_to_graph[it.first]->id;
      }
    }
  }

  void Graph::PrintDot(const std::string &path) const
  {
    if (path.empty())
    {
      return;
    }
    std::ofstream os(path);

    std::vector<std::int32_t> consensus_rank(nodes_.size(), -1);
    std::int32_t rank = 0;
    for (const auto &it : consensus_)
    {
      consensus_rank[it->id] = rank++;
    }

    os << "digraph " << sequences_.size() << " {" << std::endl
       << "  graph [rankdir = LR]" << std::endl;
    for (const auto &it : nodes_)
    {
      os << "  " << it->id << "[label = \"" << it->id << " - "
         << static_cast<char>(decoder_[it->code]) << "\"";
      if (consensus_rank[it->id] != -1)
      {
        os << ", style = filled, fillcolor = goldenrod1";
      }
      os << "]" << std::endl;

      for (const auto &jt : it->outedges)
      {
        os << "  " << it->id << " -> " << jt->head->id
           << " [label = \"" << jt->weight << "\"";
        if (consensus_rank[it->id] + 1 == consensus_rank[jt->head->id])
        {
          os << ", color = goldenrod1";
        }
        os << "]" << std::endl;
      }
      for (const auto &jt : it->aligned_nodes)
      {
        if (jt->id > it->id)
        {
          os << "  " << it->id << " -> " << jt->id
             << " [style = dotted, arrowhead = none]" << std::endl;
        }
      }
    }
    os << "}" << std::endl;

    os.close();
  }

  void Graph::Clear()
  {
    num_codes_ = 0;
    std::fill(coder_.begin(), coder_.end(), -1);
    std::fill(decoder_.begin(), decoder_.end(), -1);
    sequences_.clear();
    nodes_.clear();
    edges_.clear();
    rank_to_node_.clear();
    consensus_.clear();
  }

  // ADD by Xiao Luo
  void Graph::PruneGraph(int64_t min_weight, double min_confidence, double min_support)
  {
    // edge: u -> v
    // confidence(u->v): out-edges of u
    // confidence(v->u): in-edges  of v
    std::vector<uint64_t> edge_indexes_to_prune;
    uint64_t edge_index = 0;
    double confidence_uv;
    double confidence_vu;
    double support;
    int64_t total_weight;

    // int k = 0;
    for (const auto &it : edges_)
    {
      if (it->weight < min_weight)
      {
        // std::cerr<<"k="<<k<<"it->weight:"<<it->weight<<std::endl;
        // k++;
        edge_indexes_to_prune.emplace_back(1);
        continue;
      }
      // std::cerr<<"keeping, it->weight:"<<it->weight<<std::endl;

      //compute confidence_uv
      total_weight = 0;
      for (const auto &jt : it->tail->outedges)
      {
        total_weight += jt->weight;
      }
      confidence_uv = double(it->weight) / total_weight;

      //compute confidence_vu
      total_weight = 0;
      for (const auto &jt : it->head->inedges)
      {
        total_weight += jt->weight;
      }
      confidence_vu = double(it->weight) / total_weight;

      if (confidence_uv >= min_confidence && confidence_vu >= min_confidence && support >= min_support)
      {
        edge_indexes_to_prune.emplace_back(0);
      }
      else
      {
        edge_indexes_to_prune.emplace_back(1);
      }
    }

    std::vector<std::unique_ptr<Edge> > remained_edges;
    // int i = 0;
    for (auto &it : edges_)
    {

      if (!edge_indexes_to_prune[edge_index])
      {
        // std::cerr << "keep edge, i:" << i << std::endl;
        // std::unique_ptr<Edge> ptr = std::move(it);
        // remained_edges.emplace_back(ptr.get()); // these two lines cause a bug
        remained_edges.emplace_back(std::move(it));
      }
      else
      {
        // std::cerr << "remove edge, i:" << i << std::endl;
        for (auto &jt : it->tail->outedges)
        {
          if (jt == it.get())
          {
            jt = nullptr;
            break;
          }
        }
        for (auto &jt : it->head->inedges)
        {
          if (jt == it.get())
          {
            jt = nullptr;
            break;
          }
        }

        //remove edge
        it = nullptr;
      }
      edge_index++;
      // i++;
    }

    edges_.clear();
    edges_.shrink_to_fit();
    for (auto &it : remained_edges)
    {
      edges_.emplace_back(std::move(it));
    }
    remained_edges.clear();
    remained_edges.shrink_to_fit();
  }

  std::vector<std::uint32_t> Graph::DfsUtil(
      std::uint32_t v, bool visited[],
      std::unordered_map<Graph::Node *,
                         std::uint32_t> &nodep2idx)
  {
    // return a connected component(index of nodes) via DFS
    visited[v] = true; //v is the index in nodes_
    connected_component_.emplace_back(v);

    // recurse for all adjacent vertices of v
    std::vector<std::uint32_t> adj_vertices;
    for (auto &it : nodes_[v]->inedges)
    {
      if (it)
      {
        adj_vertices.emplace_back(nodep2idx[it->tail]); //predecessor
      }
    }

    for (auto &it : nodes_[v]->outedges)
    {
      if (it)
      {
        adj_vertices.emplace_back(nodep2idx[it->head]); //successor
      }
    }

    for (auto &u : adj_vertices)
    {
      if (!visited[u])
      {
        DfsUtil(u, visited, nodep2idx);
      }
    }
    return connected_component_;
  }

  Graph Graph::LargestSubgraph()
  {
    // extract the largest connected component as a new POA subgraph
    std::uint32_t num_nodes = nodes_.size();
    bool *visited = new bool[num_nodes];
    std::vector<std::uint32_t> connected_component;
    std::vector<std::uint32_t> largest_connected_component;

    std::uint32_t id = 0;
    std::uint32_t largest_cc_id = 0;
    std::uint32_t largest_cc_size = 0;

    std::unordered_map<Node *, std::uint32_t> nodep2idx;

    for (std::uint32_t v = 0; v < num_nodes; v++)
    {
      visited[v] = false;
      nodep2idx.insert(std::make_pair(nodes_[v].get(), v));
    }

    for (std::uint32_t v = 0; v < num_nodes; v++)
    {
      if (!visited[v])
      {
        connected_component = DfsUtil(v, visited, nodep2idx);
        connected_component_.clear();
        connected_component_.shrink_to_fit();
        // std::cerr << "connected componet size:" << connected_component.size() << std::endl;
        if (connected_component.size() >= largest_cc_size)
        {
          largest_cc_id = id;
          largest_cc_size = connected_component.size();
          largest_connected_component = connected_component;
        }
        id++;
      }
    }
    // std::cerr << "largest connected component size = " << largest_cc_size << std::endl;
    delete[] visited;

    // init subgraph
    Graph subgraph{};
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;

    //add nodes and edges for the largest subgraph
    //graph v(index in graph nodes_, not id) to sugraph v(index in subgraph nodes_, id = index)
    std::unordered_map<std::uint32_t, std::uint32_t> v2subv;
    for (const auto &v : largest_connected_component)
    {
      // std::cerr << "this is the node:" << v << std::endl;
      //
      v2subv.insert(std::make_pair(v, subgraph.AddNodeForSubgraph(nodes_[v]->code)));
    }

    for (auto &v : largest_connected_component) //O(E) because only consider outedges
    {
      for (const auto &jt : nodes_[v]->outedges)
      {
        if (jt)
        {
          subgraph.AddEdgeForSubgraph(subgraph.nodes_[v2subv[v]].get(),
                                      subgraph.nodes_[v2subv[nodep2idx[jt->head]]].get(), 0);
        }
      }
    }

    subgraph.TopologicalSort();
    return subgraph;
  }

  void Graph::AddEdgeForSubgraph(Graph::Node *tail, Graph::Node *head, std::uint32_t weight)
  {
    // std::cerr << "sequence size= " << sequences_.size() << std::endl;
    edges_.emplace_back(new Edge(tail, head, 0, weight));
    // edges_.emplace_back(new Edge(tail, head, sequences_.size(), weight));
    tail->outedges.emplace_back(edges_.back().get());
    head->inedges.emplace_back(edges_.back().get());
  }

  std::uint32_t Graph::AddNodeForSubgraph(std::uint32_t code)
  {
    // nodes_.emplace_back(new Node(v, code));
    nodes_.emplace_back(new Node(nodes_.size(), code));
    return nodes_.size() - 1;
    // return nodes_.back().get();
  }

  void Graph::AddWeights(
      const Alignment &alignment,
      const char *sequence,
      std::uint32_t sequence_len,
      const std::vector<std::uint32_t> &weights)
  {
    //align reads to the graph and add the weight of edges without adding any node
    if (sequence_len == 0)
    {
      return;
    }
    if (sequence_len != weights.size())
    {
      throw std::invalid_argument(
          "[spoa::Graph::AddWeights] error: "
          "sequence and weights are of unequal size!");
    }

    /*  for (std::uint32_t i = 0; i < sequence_len; ++i)
    {
      if (coder_[sequence[i]] == -1)
      {
        coder_[sequence[i]] = num_codes_;
        decoder_[num_codes_++] = sequence[i];
      }
    } */

    if (alignment.empty())
    {
      std::cerr << "This read fails to be aligned to the POA, skipping it..." << std::endl;
      return;
    }

    Node *prev = nullptr;
    Node *curr = nullptr;

    for (const auto &it : alignment)
    {
      // how to handle with gap ?? (-1 in the alignment)
      // only add weight of edges which are supported by raw reads (different with my prototype version)
      if (it.first == -1 || it.second == -1)
      {
        prev = nullptr;
        continue;
      }

      curr = nodes_[it.first].get(); //it.first is the id of Node, and index of nodes_
      if (prev)
      {
        // both nodes contribute to weight
        AddEdge(prev, curr, weights[it.second - 1] + weights[it.second]);
      }
      prev = curr;
    }
  }

  std::string Graph::GenerateCorrectedSequence(const Alignment &alignment)
  {
    std::string dst{};
    for (const auto &it : alignment)
    {
      if (it.first == -1)
      {
        continue;
      }
      dst += decoder_[nodes_[it.first]->code];
    }
    return dst;
  }

  // END

} // namespace spoa
