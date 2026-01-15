// graphSSR - Graph-based SSR finder
// Author: https://github.com/zergger (kai@bequ.net)
// Version: 1.0

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <thread>
#include <atomic>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kAdaptiveQuantile = 0.95;

struct SSRFragment {
    int start = 0;
    int end = 0;
    std::string motif;
    double avg_cov = 0.0;
};

struct SSRRecord {
    std::string sequence_id;
    int start = 0;
    int end = 0;
    std::string motif;
    std::string sequence;
    int size = 0;
    double avg_cov = 0.0;
};

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return base;
    }
}

std::string reverse_complement(const std::string& motif) {
    std::string rc;
    rc.reserve(motif.size());
    for (auto it = motif.rbegin(); it != motif.rend(); ++it) {
        rc.push_back(complement(*it));
    }
    return rc;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    const size_t mid = n / 2;
    if (n % 2 == 1) {
        return values[mid];
    }
    return (values[mid - 1] + values[mid]) / 2.0;
}

std::optional<double> quantile(std::vector<double> values, double q) {
    if (values.empty()) {
        return std::nullopt;
    }
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n == 1) {
        return values[0];
    }
    const double pos = q * (n - 1);
    const auto lower = static_cast<size_t>(std::floor(pos));
    const auto upper = static_cast<size_t>(std::ceil(pos));
    if (lower == upper) {
        return values[lower];
    }
    const double weight = pos - lower;
    return values[lower] + (values[upper] - values[lower]) * weight;
}

class DeBruijnGraph {
public:
    DeBruijnGraph(const std::string& sequence, int k, int min_kmer_count = 0, int max_cycle_len = -1)
        : k_(k),
          sequence_(sequence),
          min_kmer_count_(min_kmer_count),
          max_cycle_len_(max_cycle_len) {
        build();
        branch_ratio_threshold_ = auto_branch_ratio_threshold();
    }

    std::vector<SSRFragment> find_ssr_fragments(double min_avg_coverage, int max_merge_dist) {
        std::vector<SSRFragment> fragments;
        std::unordered_set<std::string> processed_cycles;

        std::vector<std::pair<std::string, int>> sorted_nodes;
        sorted_nodes.reserve(nodes_.size());
        for (const auto& kv : nodes_) {
            sorted_nodes.push_back(kv);
        }
        std::sort(sorted_nodes.begin(), sorted_nodes.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& entry : sorted_nodes) {
            const std::string& start_node = entry.first;
            const int coverage = entry.second;
            if (coverage < min_avg_coverage) {
                break;
            }

            auto cycle = find_cycle(start_node);
            if (cycle.empty()) {
                continue;
            }

            std::string motif_seq = motif_from_cycle(cycle);
            if (motif_seq.empty()) {
                continue;
            }

            std::string canonical = canonical_motif(motif_seq);
            if (canonical.empty() || processed_cycles.count(canonical) > 0) {
                continue;
            }
            processed_cycles.insert(canonical);

            double avg_cov = 0.0;
            for (const auto& node : cycle) {
                auto it = nodes_.find(node);
                if (it != nodes_.end()) {
                    avg_cov += it->second;
                }
            }
            avg_cov /= static_cast<double>(cycle.size());
            if (avg_cov < min_avg_coverage) {
                continue;
            }

            std::vector<const std::vector<int>*> position_lists;
            for (const auto& node : cycle) {
                auto it = kmer_positions_.find(node);
                if (it != kmer_positions_.end() && !it->second.empty()) {
                    position_lists.push_back(&it->second);
                }
            }
            if (position_lists.empty()) {
                continue;
            }

            auto merged_positions = merge_positions(position_lists);
            if (merged_positions.empty()) {
                continue;
            }
            const int gap = auto_cluster_gap(merged_positions, static_cast<int>(motif_seq.size()), max_merge_dist);

            std::vector<std::pair<int, int>> clusters;
            int current_start = -1;
            int current_last = -1;
            int last_pos = -1;
            for (int pos : merged_positions) {
                if (last_pos != -1 && pos == last_pos) {
                    continue;
                }
                if (current_start == -1) {
                    current_start = pos;
                    current_last = pos;
                } else if (pos - last_pos > gap) {
                    clusters.emplace_back(current_start, current_last);
                    current_start = pos;
                    current_last = pos;
                } else {
                    current_last = pos;
                }
                last_pos = pos;
            }
            if (current_start != -1) {
                clusters.emplace_back(current_start, current_last);
            }

            for (const auto& cluster : clusters) {
                int end_pos = cluster.second + k_;
                fragments.push_back({cluster.first + 1, end_pos, canonical, avg_cov});
            }
        }

        return fragments;
    }

private:
    int k_;
    std::string sequence_;
    int min_kmer_count_;
    int max_cycle_len_;
    std::unordered_map<std::string, int> nodes_;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> graph_;
    std::unordered_map<std::string, std::vector<int>> kmer_positions_;
    std::optional<double> branch_ratio_threshold_;

    void build() {
        if (static_cast<int>(sequence_.size()) < k_) {
            return;
        }

        if (min_kmer_count_ <= 0) {
            std::unordered_map<std::string, int> counts;
            for (int i = 0; i <= static_cast<int>(sequence_.size()) - k_; ++i) {
                std::string kmer = sequence_.substr(i, k_);
                if (kmer.find('N') != std::string::npos) {
                    continue;
                }
                counts[kmer] += 1;
            }
            nodes_ = std::move(counts);
            min_kmer_count_ = auto_min_kmer_count(nodes_);
            for (int i = 0; i <= static_cast<int>(sequence_.size()) - k_; ++i) {
                std::string kmer = sequence_.substr(i, k_);
                if (kmer.find('N') != std::string::npos) {
                    continue;
                }
                auto it = nodes_.find(kmer);
                if (it == nodes_.end() || it->second < min_kmer_count_) {
                    continue;
                }
                kmer_positions_[kmer].push_back(i);
                if (i < static_cast<int>(sequence_.size()) - k_) {
                    std::string next_kmer = sequence_.substr(i + 1, k_);
                    if (next_kmer.find('N') == std::string::npos) {
                        auto it_next = nodes_.find(next_kmer);
                        if (it_next != nodes_.end() && it_next->second >= min_kmer_count_) {
                            graph_[kmer][next_kmer] += 1;
                        }
                    }
                }
            }
        } else if (min_kmer_count_ > 1) {
            for (int i = 0; i <= static_cast<int>(sequence_.size()) - k_; ++i) {
                std::string kmer = sequence_.substr(i, k_);
                if (kmer.find('N') != std::string::npos) {
                    continue;
                }
                nodes_[kmer] += 1;
            }
            for (int i = 0; i <= static_cast<int>(sequence_.size()) - k_; ++i) {
                std::string kmer = sequence_.substr(i, k_);
                if (kmer.find('N') != std::string::npos) {
                    continue;
                }
                if (nodes_[kmer] < min_kmer_count_) {
                    continue;
                }
                kmer_positions_[kmer].push_back(i);
                if (i < static_cast<int>(sequence_.size()) - k_) {
                    std::string next_kmer = sequence_.substr(i + 1, k_);
                    if (next_kmer.find('N') == std::string::npos && nodes_[next_kmer] >= min_kmer_count_) {
                        graph_[kmer][next_kmer] += 1;
                    }
                }
            }
        } else {
            for (int i = 0; i <= static_cast<int>(sequence_.size()) - k_; ++i) {
                std::string kmer = sequence_.substr(i, k_);
                if (kmer.find('N') != std::string::npos) {
                    continue;
                }
                nodes_[kmer] += 1;
                kmer_positions_[kmer].push_back(i);
                if (i < static_cast<int>(sequence_.size()) - k_) {
                    std::string next_kmer = sequence_.substr(i + 1, k_);
                    if (next_kmer.find('N') == std::string::npos) {
                        graph_[kmer][next_kmer] += 1;
                    }
                }
            }
        }
    }

    std::optional<double> auto_branch_ratio_threshold() const {
        std::vector<double> ratios;
        for (const auto& kv : graph_) {
            if (kv.second.size() < 3) {
                continue;
            }
            std::vector<int> weights;
            weights.reserve(kv.second.size());
            for (const auto& edge : kv.second) {
                weights.push_back(edge.second);
            }
            std::sort(weights.begin(), weights.end(), std::greater<int>());
            const int second_weight = weights[1];
            const int third_weight = weights[2];
            if (second_weight <= 0) {
                continue;
            }
            ratios.push_back(static_cast<double>(third_weight) / static_cast<double>(second_weight));
        }
        if (ratios.empty()) {
            return std::nullopt;
        }
        const double med = median(ratios);
        std::vector<double> deviations;
        deviations.reserve(ratios.size());
        for (double r : ratios) {
            deviations.push_back(std::abs(r - med));
        }
        const double mad = median(deviations);
        double threshold;
        if (mad > 0) {
            threshold = med + 3 * mad;
        } else {
            auto q = quantile(ratios, kAdaptiveQuantile);
            if (!q.has_value()) {
                return std::nullopt;
            }
            threshold = q.value();
        }
        threshold = std::max(0.0, std::min(1.0, threshold));
        return threshold;
    }

    int auto_min_kmer_count(const std::unordered_map<std::string, int>& counts) const {
        if (counts.empty()) {
            return 1;
        }
        std::vector<double> values;
        values.reserve(counts.size());
        for (const auto& kv : counts) {
            values.push_back(static_cast<double>(kv.second));
        }
        const double med = median(values);
        std::vector<double> deviations;
        deviations.reserve(values.size());
        for (double v : values) {
            deviations.push_back(std::abs(v - med));
        }
        const double mad = median(deviations);
        double threshold;
        if (mad > 0) {
            threshold = med + 3 * mad;
        } else {
            auto q = quantile(values, kAdaptiveQuantile);
            threshold = q.value_or(1.0);
        }
        return std::max(1, static_cast<int>(std::ceil(threshold)));
    }

    std::vector<int> merge_positions(const std::vector<const std::vector<int>*>& lists) const {
        struct HeapEntry {
            int pos;
            size_t list_idx;
            size_t elem_idx;
        };
        struct Compare {
            bool operator()(const HeapEntry& a, const HeapEntry& b) const {
                return a.pos > b.pos;
            }
        };
        std::priority_queue<HeapEntry, std::vector<HeapEntry>, Compare> heap;
        for (size_t i = 0; i < lists.size(); ++i) {
            if (!lists[i]->empty()) {
                heap.push({(*lists[i])[0], i, 0});
            }
        }
        std::vector<int> merged;
        merged.reserve(128);
        int last_pos = std::numeric_limits<int>::min();
        while (!heap.empty()) {
            auto current = heap.top();
            heap.pop();
            if (current.pos != last_pos) {
                merged.push_back(current.pos);
                last_pos = current.pos;
            }
            const auto& list = *lists[current.list_idx];
            const size_t next_idx = current.elem_idx + 1;
            if (next_idx < list.size()) {
                heap.push({list[next_idx], current.list_idx, next_idx});
            }
        }
        return merged;
    }

    int auto_cluster_gap(const std::vector<int>& merged_positions, int motif_len, int max_merge_dist) const {
        if (merged_positions.empty()) {
            return motif_len;
        }
        std::vector<double> gaps;
        gaps.reserve(merged_positions.size());
        int last_pos = merged_positions.front();
        for (size_t i = 1; i < merged_positions.size(); ++i) {
            int gap = merged_positions[i] - last_pos;
            if (gap > 0) {
                gaps.push_back(static_cast<double>(gap));
            }
            last_pos = merged_positions[i];
        }
        int gap = motif_len;
        if (!gaps.empty()) {
            const double med = median(gaps);
            std::vector<double> deviations;
            deviations.reserve(gaps.size());
            for (double g : gaps) {
                deviations.push_back(std::abs(g - med));
            }
            const double mad = median(deviations);
            double estimate;
            if (mad > 0) {
                estimate = med + 3 * mad;
            } else {
                auto q = quantile(gaps, kAdaptiveQuantile);
                estimate = q.value_or(static_cast<double>(motif_len));
            }
            gap = static_cast<int>(std::ceil(estimate));
        }
        gap = std::max(1, gap);
        if (max_merge_dist > 0) {
            gap = std::min(gap, max_merge_dist);
        }
        return gap;
    }

    int adaptive_cycle_len(const std::string& start_node) const {
        auto it = kmer_positions_.find(start_node);
        if (it == kmer_positions_.end() || it->second.size() < 3) {
            return k_ + 7;
        }
        const auto& positions = it->second;
        std::vector<double> gaps;
        gaps.reserve(positions.size());
        int last_pos = positions[0];
        for (size_t i = 1; i < positions.size(); ++i) {
            int gap = positions[i] - last_pos;
            if (gap > 0) {
                gaps.push_back(static_cast<double>(gap));
            }
            last_pos = positions[i];
        }
        if (gaps.size() < 2) {
            return k_ + 7;
        }
        const double med = median(gaps);
        std::vector<double> deviations;
        deviations.reserve(gaps.size());
        for (double g : gaps) {
            deviations.push_back(std::abs(g - med));
        }
        const double mad = median(deviations);
        double estimate;
        if (mad > 0) {
            estimate = med + 3 * mad;
        } else {
            auto q = quantile(gaps, kAdaptiveQuantile);
            estimate = q.value_or(static_cast<double>(k_ + 7));
        }
        if (estimate <= 0) {
            return k_ + 7;
        }
        return std::max(1, static_cast<int>(std::ceil(estimate)));
    }

    std::vector<std::pair<std::string, int>> select_successors(
        const std::unordered_map<std::string, int>& successors,
        int max_branches,
        std::optional<double> adaptive_ratio) const {
        std::vector<std::pair<std::string, int>> ranked;
        ranked.reserve(successors.size());
        for (const auto& kv : successors) {
            ranked.push_back(kv);
        }
        std::sort(ranked.begin(), ranked.end(),
                  [](const auto& a, const auto& b) {
                      if (a.second != b.second) {
                          return a.second > b.second;
                      }
                      return a.first < b.first;
                  });
        int limit = max_branches;
        if (adaptive_ratio.has_value() && ranked.size() >= 3) {
            const int second_weight = ranked[1].second;
            const int third_weight = ranked[2].second;
            if (second_weight > 0 && third_weight >= adaptive_ratio.value() * second_weight) {
                limit = std::max(limit, 3);
            }
        }
        if (static_cast<int>(ranked.size()) > limit) {
            ranked.resize(static_cast<size_t>(limit));
        }
        return ranked;
    }

    std::vector<std::string> find_cycle(const std::string& start_node) {
        int max_len = max_cycle_len_;
        if (max_len <= 0) {
            max_len = adaptive_cycle_len(start_node);
        }
        auto cycle = find_cycle_with_branches(start_node, max_len, 2, std::nullopt);
        if (!cycle.empty()) {
            return cycle;
        }
        if (!branch_ratio_threshold_.has_value()) {
            return {};
        }
        return find_cycle_with_branches(start_node, max_len, 2, branch_ratio_threshold_);
    }

    std::vector<std::string> find_cycle_with_branches(
        const std::string& start_node,
        int max_len,
        int max_branches,
        std::optional<double> adaptive_ratio) const {
        struct StackNode {
            std::string current;
            std::vector<std::string> path;
            std::unordered_set<std::string> visited;
        };
        std::vector<StackNode> stack;
        stack.push_back({start_node, {start_node}, {start_node}});
        while (!stack.empty()) {
            auto node = std::move(stack.back());
            stack.pop_back();
            if (static_cast<int>(node.path.size()) > max_len) {
                continue;
            }
            auto it = graph_.find(node.current);
            if (it == graph_.end() || it->second.empty()) {
                continue;
            }
            auto ranked = select_successors(it->second, max_branches, adaptive_ratio);
            for (const auto& next_pair : ranked) {
                const std::string& next_node = next_pair.first;
                if (node.visited.count(next_node)) {
                    auto pos = std::find(node.path.begin(), node.path.end(), next_node);
                    if (pos != node.path.end()) {
                        return std::vector<std::string>(pos, node.path.end());
                    }
                    continue;
                }
                if (static_cast<int>(node.path.size()) < max_len) {
                    StackNode next;
                    next.current = next_node;
                    next.path = node.path;
                    next.path.push_back(next_node);
                    next.visited = node.visited;
                    next.visited.insert(next_node);
                    stack.push_back(std::move(next));
                }
            }
        }
        return {};
    }

    std::string motif_from_cycle(const std::vector<std::string>& cycle) const {
        std::string motif;
        motif.reserve(cycle.size());
        for (const auto& node : cycle) {
            if (!node.empty()) {
                motif.push_back(node[0]);
            }
        }
        return motif;
    }

    std::string canonical_motif(const std::string& motif) const {
        if (motif.empty()) {
            return "";
        }
        std::string rc = reverse_complement(motif);
        std::string best;
        for (const auto& candidate : {motif, rc}) {
            for (size_t i = 0; i < candidate.size(); ++i) {
                std::string rotated = candidate.substr(i) + candidate.substr(0, i);
                if (best.empty() || rotated < best) {
                    best = std::move(rotated);
                }
            }
        }
        return best;
    }
};

std::vector<SSRRecord> find_ssrs_in_sequence(
    const std::string& sequence,
    const std::string& sequence_id,
    int k,
    int min_len,
    double min_cov,
    int max_merge_dist,
    bool verbose) {
    if (verbose) {
        std::cout << "\n--- Processing sequence: " << sequence_id
                  << " (length: " << sequence.size() << " bp) ---\n";
        std::cout << "  - Step 1: building de Bruijn graph (k=" << k << ")...\n";
    }

    DeBruijnGraph graph(sequence, k);

    if (verbose) {
        std::cout << "  - Step 2: finding SSR fragments...\n";
    }
    auto ssr_fragments = graph.find_ssr_fragments(min_cov, max_merge_dist);

    if (verbose) {
        std::cout << "  - Step 3: merging and trimming " << ssr_fragments.size() << " fragments...\n";
    }
    if (ssr_fragments.empty()) {
        if (verbose) {
            std::cout << "  - Step 4: no SSR found.\n";
        }
        return {};
    }

    std::sort(ssr_fragments.begin(), ssr_fragments.end(),
              [](const SSRFragment& a, const SSRFragment& b) { return a.start < b.start; });

    std::vector<SSRFragment> merged;
    SSRFragment current = ssr_fragments.front();
    for (size_t i = 1; i < ssr_fragments.size(); ++i) {
        const auto& next = ssr_fragments[i];
        const int dist = next.start - current.end;
        if (next.motif == current.motif && dist <= max_merge_dist) {
            current.end = std::max(current.end, next.end);
        } else {
            merged.push_back(current);
            current = next;
        }
    }
    merged.push_back(current);

    std::vector<SSRRecord> final_ssrs;
    for (auto& ssr : merged) {
        int size = ssr.end - ssr.start + 1;
        if (size < min_len) {
            continue;
        }
        const int motif_len = static_cast<int>(ssr.motif.size());
        if (motif_len == 0) {
            continue;
        }
        int new_size = (size / motif_len) * motif_len;
        if (new_size < motif_len) {
            continue;
        }
        const int new_end = ssr.start + new_size - 1;
        if (new_end > static_cast<int>(sequence.size())) {
            continue;
        }
        SSRRecord record;
        record.sequence_id = sequence_id;
        record.start = ssr.start;
        record.end = new_end;
        record.motif = ssr.motif;
        record.size = new_size;
        record.avg_cov = ssr.avg_cov;
        record.sequence = sequence.substr(static_cast<size_t>(record.start - 1),
                                          static_cast<size_t>(record.size));
        final_ssrs.push_back(std::move(record));
    }

    if (verbose) {
        std::cout << "  - Step 4: found " << final_ssrs.size() << " SSRs.\n";
    }
    return final_ssrs;
}

bool read_fasta(const std::string& path,
                std::vector<std::pair<std::string, std::string>>& records) {
    std::ifstream input(path);
    if (!input) {
        return false;
    }
    std::string line;
    std::string current_id;
    std::ostringstream current_seq;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!current_id.empty()) {
                records.emplace_back(current_id, current_seq.str());
                current_seq.str(std::string());
                current_seq.clear();
            }
            std::istringstream iss(line.substr(1));
            iss >> current_id;
        } else {
            current_seq << line;
        }
    }
    if (!current_id.empty()) {
        records.emplace_back(current_id, current_seq.str());
    }
    return true;
}

void write_results(const std::vector<SSRRecord>& records, const std::string& output_filename) {
    if (records.empty()) {
        std::cout << "\nNo SSRs found; output file not generated.\n";
        return;
    }
    std::cout << "\n--- Writing output file ---\n";

    std::vector<SSRRecord> sorted_records = records;
    std::sort(sorted_records.begin(), sorted_records.end(),
              [](const SSRRecord& a, const SSRRecord& b) {
                  if (a.sequence_id != b.sequence_id) {
                      return a.sequence_id < b.sequence_id;
                  }
                  return a.start < b.start;
              });

    std::ofstream out(output_filename);
    out << "Sequence_ID\tSSR_nr\tInferred_Motif\tSSR_Sequence\tSize\tStart\tEnd\tAvg_Coverage\n";
    int idx = 1;
    for (const auto& ssr : sorted_records) {
        out << ssr.sequence_id << '\t'
            << idx++ << '\t'
            << ssr.motif << '\t'
            << ssr.sequence << '\t'
            << ssr.size << '\t'
            << ssr.start << '\t'
            << ssr.end << '\t'
            << std::fixed << std::setprecision(2) << ssr.avg_cov
            << '\n';
    }
    std::cout << "  - Output saved to: " << output_filename << "\n";
}

void print_usage(const char* prog) {
    std::cerr << "graphSSR v1.0\n"
              << "Author: https://github.com/zergger (kai@bequ.net)\n"
              << "Usage: " << prog
              << " <input_fasta> <output_prefix> [-k <kmer>] "
              << "[--min_len <len>] [--min_cov <cov>] [--max_merge_dist <dist>] "
              << "[--threads <n>] [--quiet]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    std::string input_fasta = argv[1];
    std::string output_prefix = argv[2];
    int kmer_size = 21;
    int min_len = 30;
    double min_cov = 2.0;
    int max_merge_dist = 25;
    int threads = 1;
    bool quiet = false;

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-k" || arg == "--kmer_size") && i + 1 < argc) {
            kmer_size = std::stoi(argv[++i]);
        } else if (arg == "--min_len" && i + 1 < argc) {
            min_len = std::stoi(argv[++i]);
        } else if (arg == "--min_cov" && i + 1 < argc) {
            min_cov = std::stod(argv[++i]);
        } else if (arg == "--max_merge_dist" && i + 1 < argc) {
            max_merge_dist = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (threads == 0) {
        threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    }
    if (threads < 1) {
        std::cerr << "Error: --threads must be >= 1 or 0 (auto).\n";
        return 1;
    }
    const bool verbose = !quiet && threads == 1;

    auto start = std::chrono::steady_clock::now();
    std::cout << "\n==================================================\n";
    std::cout << "Analyzing file: " << input_fasta << "\n";
    if (threads > 1) {
        std::cout << "Parallel processing with threads=" << threads << "\n";
    }
    std::cout << "==================================================\n";

    std::vector<std::pair<std::string, std::string>> sequences;
    if (!read_fasta(input_fasta, sequences)) {
        std::cerr << "Error: input file not found -> " << input_fasta << "\n";
        return 1;
    }

    std::vector<SSRRecord> all_ssrs;
    if (threads == 1) {
        for (const auto& record : sequences) {
            const auto& seq_id = record.first;
            const auto& seq = record.second;
            auto ssrs = find_ssrs_in_sequence(seq, seq_id, kmer_size, min_len, min_cov, max_merge_dist, verbose);
            all_ssrs.insert(all_ssrs.end(), ssrs.begin(), ssrs.end());
        }
    } else {
        std::vector<std::vector<SSRRecord>> per_seq_results(sequences.size());
        std::atomic<size_t> next_idx{0};
        auto worker = [&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= sequences.size()) {
                    break;
                }
                const auto& seq_id = sequences[idx].first;
                const auto& seq = sequences[idx].second;
                per_seq_results[idx] = find_ssrs_in_sequence(
                    seq, seq_id, kmer_size, min_len, min_cov, max_merge_dist, false
                );
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(threads));
        for (int t = 0; t < threads; ++t) {
            pool.emplace_back(worker);
        }
        for (auto& th : pool) {
            th.join();
        }
        for (auto& vec : per_seq_results) {
            all_ssrs.insert(all_ssrs.end(), vec.begin(), vec.end());
        }
    }

    std::string output_filename = output_prefix + ".ssrs";
    write_results(all_ssrs, output_filename);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\n==================================================\n";
    std::cout << "Done. Total sequences: " << sequences.size()
              << ", total SSRs: " << all_ssrs.size() << "\n";
    std::cout << "Output: " << output_filename << "\n";
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    std::cout << "==================================================\n";
    return 0;
}
