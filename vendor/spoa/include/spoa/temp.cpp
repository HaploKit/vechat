
std::vector<std::unique_ptr<Graph::Node> > Graph::DfsUtil(std::unique_ptr<Graph::Node> &v, bool visited[])
// std::vector<const Graph::Node *> Graph::DfsUtil(const Graph::Node *v, bool visited[])
{
    // std::unique_ptr<Graph::Node> v(w.release());
    std::vector<std::unique_ptr<Node> > connected_component;
    // std::vector<const Node *> connected_component;
    visited[v->id] = true;

    // recurse for all adjacent vertices of v
    std::vector<std::unique_ptr<Node> > adj_vertices;
    // std::vector<const Node *> adj_vertices;

    std::cerr << "\n#### v->id: " << v->id << std::endl;
    for (auto &it : v->inedges)
    {
        if (it)
        {
            std::cerr << "\ninedge(-->):" << it->tail->id << "\t" << it->head->id << std::endl;
            adj_vertices.emplace_back(std::move(it->tail)); //predecessor
        }
    }
    for (auto &it : v->outedges)
    {
        if (it)
        {
            std::cerr << it->tail << std::endl;
            std::cerr << "\noutedge(-->):" << it->tail->id << "\t" << it->head->id << std::endl;
            adj_vertices.emplace_back(std::move(it->head)); //successor
        }
        else
        {
            std::cerr << "XXXXXXXXXXXXXXX\n";
        }
    }

    std::cerr << "\nOutput node:  " << v->id << std::endl;
    connected_component.emplace_back(std::move(v));
    std::cerr << "\nOutput size:  " << connected_component.size() << std::endl;
    std::cerr << "\n #adj nodes  " << adj_vertices.size() << std::endl;

    if (!adj_vertices.empty())
    {
        for (auto &u : adj_vertices)
        {
            std::cerr << " adj node id: " << u->id;
            if (!visited[u->id])
            {
                std::cerr << " \nlast xxxx: " << u->id;
                DfsUtil(u, visited);
            }
        }
    }
    std::cerr << "\n####break###\n";
    return connected_component;
}

Graph Graph::LargestSubgraph()
{
    // extract the largest connected component as a new POA graph
    std::uint32_t num_nodes = nodes_.size();
    bool *visited = new bool[num_nodes];
    // std::vector<const Node *> connected_component;
    // std::vector<const Node *> largest_connected_component;
    // std::vector<std::vector<const Node *>> connected_components;
    std::vector<std::unique_ptr<Node> > connected_component;
    // std::vector< std::unique_ptr<Node>> largest_connected_component;
    std::vector<std::vector<std::unique_ptr<Node> > *> connected_components;

    std::uint32_t id = 0;
    std::uint32_t largest_cc_id = 0;
    std::uint32_t largest_cc_size = 0;

    for (std::uint32_t i = 0; i < num_nodes; i++)
        visited[i] = false;
    for (auto &v : nodes_)
    {
        std::cerr << " running..." << std::endl;
        if (!visited[v->id])
        {
            std::cerr << " go into if..." << std::endl;
            DfsUtil(v, visited);
            // connected_component = DfsUtil(v, visited);
            std::cerr << "connected_component running..." << std::endl;

            connected_components.emplace_back(std::move(&connected_component));
            if (connected_component.size() > largest_cc_size)
            {
                largest_cc_id = id;
                largest_cc_size = connected_component.size();
            }
            id++;
        }
    }
    // largest_connected_component = *connected_components[largest_cc_id];

    // init subgraph
    Graph subgraph{};
    subgraph.num_codes_ = num_codes_;
    subgraph.coder_ = coder_;
    subgraph.decoder_ = decoder_;

    //traverse each node in the largest connected component (subgraph)
    for (const auto &it : (*connected_components[largest_cc_id]))
    {
        subgraph.AddNode(it->code);
        for (const auto &jt : it->inedges)
        {
            subgraph.AddEdge(jt->tail, jt->head, 0);
        }
        for (const auto &jt : it->outedges)
        {
            subgraph.AddEdge(jt->tail, jt->head, 0);
        }
    }

    subgraph.TopologicalSort();

    return subgraph;
}


for (const auto &v : largest_connected_component)
    {
      std::cerr << "this is the node:" << v << std::endl;
      subgraph.AddNode(nodes_[v]->code);
      for (const auto &jt : nodes_[v]->inedges)
      {
        if (jt )
        {
          std::cerr << "hello\n";
          subgraph.AddEdgeForSubgraph(jt->tail, jt->head, 0);
        }
      }

      for (const auto &jt : nodes_[v]->outedges)
      {
        if (jt)
        {
          std::cerr << "hello\n";
          subgraph.AddEdgeForSubgraph(jt->tail, jt->head,0);
        }
      }
    }