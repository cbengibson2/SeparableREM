// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// Build sorted adjacency lists (out-neighbours and in-neighbours) from an
// edgelist.  Nodes are represented by their raw integer IDs so no string
// conversion is needed.
//
// Returns a named list with two elements:
//   $out  – List, one IntegerVector of out-neighbours per node ID
//   $in   – List, one IntegerVector of in-neighbours per node ID
//
// [[Rcpp::export]]
List build_adjacency(IntegerVector senders, IntegerVector receivers) {
  int m = senders.size();

  // Collect unique node IDs and map them to 0-based indices
  std::vector<int> node_ids;
  std::unordered_map<int, int> id_to_idx;
  for (int i = 0; i < m; ++i) {
    int s = senders[i], r = receivers[i];
    if (id_to_idx.find(s) == id_to_idx.end()) {
      id_to_idx[s] = node_ids.size();
      node_ids.push_back(s);
    }
    if (id_to_idx.find(r) == id_to_idx.end()) {
      id_to_idx[r] = node_ids.size();
      node_ids.push_back(r);
    }
  }
  int n = node_ids.size();

  std::vector<std::vector<int>> out_adj(n), in_adj(n);
  for (int i = 0; i < m; ++i) {
    int si = id_to_idx[senders[i]];
    int ri = id_to_idx[receivers[i]];
    out_adj[si].push_back(node_ids[ri]);
    in_adj[ri].push_back(node_ids[si]);
  }

  // Deduplicate and sort each neighbour list so std::set_intersection works
  List out_list(n), in_list(n);
  CharacterVector node_names(n);
  for (int i = 0; i < n; ++i) {
    std::sort(out_adj[i].begin(), out_adj[i].end());
    out_adj[i].erase(std::unique(out_adj[i].begin(), out_adj[i].end()),
                     out_adj[i].end());
    std::sort(in_adj[i].begin(), in_adj[i].end());
    in_adj[i].erase(std::unique(in_adj[i].begin(), in_adj[i].end()),
                    in_adj[i].end());

    out_list[i] = wrap(out_adj[i]);
    in_list[i]  = wrap(in_adj[i]);
    node_names[i] = std::to_string(node_ids[i]);
  }
  out_list.attr("names") = node_names;
  in_list.attr("names")  = node_names;

  return List::create(Named("out") = out_list, Named("in") = in_list);
}

// Count shared neighbours for each (sender, receiver) pair in the risk set.
//
// adj       – output of build_adjacency(): list with $out and $in
// rs_senders, rs_receivers – integer vectors of the risk-set dyads
//
// Returns an integer matrix with columns OSP, ISP, ITP, OTP.
//
// [[Rcpp::export]]
IntegerMatrix triadic_counts(List adj,
                              IntegerVector rs_senders,
                              IntegerVector rs_receivers) {
  List out_adj = adj["out"];
  List in_adj  = adj["in"];

  // Build a fast node-ID -> neighbour-vector lookup
  CharacterVector node_names = out_adj.attr("names");
  std::unordered_map<int, std::vector<int>> out_map, in_map;
  int nn = node_names.size();
  for (int i = 0; i < nn; ++i) {
    int id = std::stoi(std::string(node_names[i]));
    IntegerVector ov = out_adj[i];
    IntegerVector iv = in_adj[i];
    out_map[id] = std::vector<int>(ov.begin(), ov.end());
    in_map[id]  = std::vector<int>(iv.begin(), iv.end());
  }

  static const std::vector<int> empty;
  auto get_out = [&](int id) -> const std::vector<int>& {
    auto it = out_map.find(id); return it != out_map.end() ? it->second : empty;
  };
  auto get_in = [&](int id) -> const std::vector<int>& {
    auto it = in_map.find(id);  return it != in_map.end()  ? it->second : empty;
  };

  int R = rs_senders.size();
  IntegerMatrix result(R, 4);
  colnames(result) = CharacterVector::create("OSP", "ISP", "ITP", "OTP");

  std::vector<int> buf;
  buf.reserve(64);

  for (int i = 0; i < R; ++i) {
    int a = rs_senders[i], b = rs_receivers[i];

    const std::vector<int>& a_out = get_out(a);
    const std::vector<int>& a_in  = get_in(a);
    const std::vector<int>& b_out = get_out(b);
    const std::vector<int>& b_in  = get_in(b);

    // OSP: shared out-neighbours of a and b
    buf.clear();
    std::set_intersection(a_out.begin(), a_out.end(),
                          b_out.begin(), b_out.end(),
                          std::back_inserter(buf));
    result(i, 0) = buf.size();

    // ISP: shared in-neighbours of a and b
    buf.clear();
    std::set_intersection(a_in.begin(), a_in.end(),
                          b_in.begin(), b_in.end(),
                          std::back_inserter(buf));
    result(i, 1) = buf.size();

    // ITP: in-neighbours of a that are out-neighbours of b
    buf.clear();
    std::set_intersection(a_in.begin(), a_in.end(),
                          b_out.begin(), b_out.end(),
                          std::back_inserter(buf));
    result(i, 2) = buf.size();

    // OTP: out-neighbours of a that are in-neighbours of b
    buf.clear();
    std::set_intersection(a_out.begin(), a_out.end(),
                          b_in.begin(), b_in.end(),
                          std::back_inserter(buf));
    result(i, 3) = buf.size();
  }

  return result;
}
