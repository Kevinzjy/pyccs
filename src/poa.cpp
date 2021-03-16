#include <string.h>
#include "spoa/spoa.hpp"

extern "C" {
    const char* poa(char** seqs, int num_seqs,
        int l, int m, int n, int g, int e, int q, int c) {

        std::vector<std::string> sequences;
        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
        }

        auto alignment_engine = spoa::AlignmentEngine::Create(static_cast<spoa::AlignmentType>(l),
            m, n, g, e, q, c);

        spoa::Graph graph{};
        
        for (const auto& it: sequences) {
            auto alignment = alignment_engine->Align(it, graph);
            graph.AddAlignment(alignment, it);
        }

        auto cns = graph.GenerateConsensus();
        
        char *consensus;
        consensus = new char [cns.size() + 1];
        strcpy (consensus, cns.c_str());
        return consensus;
    }
}
