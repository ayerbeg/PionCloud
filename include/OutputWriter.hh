#pragma once
// ============================================================
//  OutputWriter.hh
//  Pion Cloud Model -- ASCII + ROOT output
//
//  OutputWriter::write() dispatches to the sub-writers
//  according to the OutputMode bitmask in RunParams.
//
//  ROOT output is conditionally compiled when HAVE_ROOT is
//  defined at build time. If ROOT is not available but
//  OUT_ROOT is requested, a warning is printed and that
//  part is skipped gracefully.
// ============================================================

#include "PhysicsParams.hh"
#include "Results.hh"
#include <string>

namespace PionCloud {

class OutputWriter {
public:
    OutputWriter(const RunParams &p, const std::string &baseName);

    // Dispatch to ASCII and/or ROOT writers based on p_.outputMode
    void write(const ScanResults &results) const;

private:
    void writeASCII(const ScanResults &results) const;
    void writeROOT (const ScanResults &results) const;

    RunParams   p_;
    std::string base_;
};

} // namespace PionCloud
