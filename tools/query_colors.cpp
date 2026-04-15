#include <fstream>
#include <string_view>

#include "external/sshash/include/streaming_query.hpp"

using namespace fulgor;

template <typename FulgorIndex>
void print_colors(std::ostream& out, FulgorIndex const& index, uint32_t color_set_id) {
    auto it = index.color_set(color_set_id);
    out << it.size();
    for (uint64_t i = 0; i != it.size(); ++i, ++it) out << '\t' << *it;
}

template <typename FulgorIndex>
void print_lookup_result(std::ostream& out, FulgorIndex const& index, uint64_t pos,
                         std::string_view kmer, sshash::lookup_result const& answer) {
    out << pos << '\t' << kmer << '\t';
    if (answer.kmer_id == sshash::constants::invalid_uint64) {
        out << 0 << '\t' << -1 << '\t' << 0 << '\n';
        return;
    }

    uint32_t color_set_id = index.u2c(answer.contig_id);
    out << int(answer.kmer_orientation) << '\t' << color_set_id << '\t';
    print_colors(out, index, color_set_id);
    out << '\n';
}

template <typename FulgorIndex>
void query_sequence(std::ostream& out, FulgorIndex const& index, std::string const& sequence,
                    bool strand_specific) {
    const uint64_t k = index.k();
    const uint64_t num_kmers = sequence.length() - k + 1;
    auto emit = [&](uint64_t pos, sshash::lookup_result const& answer) {
        print_lookup_result(out, index, pos, std::string_view(sequence.data() + pos, k), answer);
    };

    if (index.get_k2u().canonical()) {
        sshash::streaming_query<kmer_type, true> query(&index.get_k2u());
        query.reset();
        for (uint64_t i = 0; i != num_kmers; ++i) emit(i, query.lookup_advanced(sequence.data() + i));
    } else if (strand_specific) {
        sshash::streaming_query<kmer_type, false> query(&index.get_k2u(), false);
        query.reset();
        for (uint64_t i = 0; i != num_kmers; ++i) emit(i, query.lookup_advanced(sequence.data() + i));
    } else {
        sshash::streaming_query<kmer_type, false> query(&index.get_k2u());
        query.reset();
        for (uint64_t i = 0; i != num_kmers; ++i) emit(i, query.lookup_advanced(sequence.data() + i));
    }
}

template <typename FulgorIndex>
int query_colors(std::string const& index_filename, std::string const& output_filename,
                 std::string const& kmer, std::string const& sequence, bool strand_specific,
                 bool verbose) {
    FulgorIndex index;
    if (verbose) essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    if (verbose) essentials::logger("DONE");

    if (strand_specific && index.get_k2u().canonical()) {
        std::cerr << "Error: --strand-specific requires a non-canonical hybrid index."
                  << std::endl;
        return 1;
    }

    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!output_filename.empty()) {
        out_file.open(output_filename, std::ios::out | std::ios::trunc);
        if (!out_file) {
            std::cerr << "could not open output file " + output_filename << std::endl;
            return 1;
        }
        out = &out_file;
    }

    if (!kmer.empty()) {
        if (kmer.length() != index.k()) {
            std::cerr << "Error: --kmer must have length exactly k = " << index.k() << '.'
                      << std::endl;
            return 1;
        }
        sshash::lookup_result answer;
        if (sshash::util::is_valid<kmer_type>(kmer.c_str(), kmer.length())) {
            answer = index.get_k2u().lookup_advanced(kmer.c_str(), !strand_specific);
        }
        print_lookup_result(*out, index, 0, kmer, answer);
        return 0;
    }

    if (sequence.length() < index.k()) {
        std::cerr << "Error: --sequence must have length at least k = " << index.k() << '.'
                  << std::endl;
        return 1;
    }

    query_sequence(*out, index, sequence, strand_specific);
    return 0;
}

int query_colors(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("kmer", "Query one k-mer and print its full color set.", "--kmer", false);
    parser.add("sequence", "Query every k-mer in a sequence and print all color sets.", "--sequence",
               false);
    parser.add("output_filename",
               "File where output will be written. Default is stdout when -o is omitted.", "-o",
               false);
    parser.add("strand_specific",
               "Disable reverse-complement fallback. Requires a non-canonical hybrid index.",
               "--strand-specific", false, true);
    parser.add("verbose", "Verbose output during query (default is false).", "--verbose", false,
               true);
    if (!parser.parse()) return 1;

    bool has_kmer = parser.parsed("kmer");
    bool has_sequence = parser.parsed("sequence");
    if (has_kmer == has_sequence) {
        std::cerr << "Specify exactly one of --kmer or --sequence." << std::endl;
        return 1;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    auto output_filename =
        parser.parsed("output_filename") ? parser.get<std::string>("output_filename") : "";
    auto kmer = has_kmer ? parser.get<std::string>("kmer") : "";
    auto sequence = has_sequence ? parser.get<std::string>("sequence") : "";
    bool strand_specific = parser.get<bool>("strand_specific");

    if (is_meta(index_filename)) {
        return query_colors<mfur_index_t>(index_filename, output_filename, kmer, sequence,
                                          strand_specific, verbose);
    } else if (is_meta_diff(index_filename)) {
        return query_colors<mdfur_index_t>(index_filename, output_filename, kmer, sequence,
                                           strand_specific, verbose);
    } else if (is_diff(index_filename)) {
        return query_colors<dfur_index_t>(index_filename, output_filename, kmer, sequence,
                                          strand_specific, verbose);
    } else if (is_hybrid(index_filename)) {
        return query_colors<hfur_index_t>(index_filename, output_filename, kmer, sequence,
                                          strand_specific, verbose);
    }

    std::cerr << "Wrong index filename supplied." << std::endl;
    return 1;
}
