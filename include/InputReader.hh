#pragma once
// ============================================================
//  InputReader.hh
//  Pion Cloud Model -- input file parser
//
//  PURPOSE
//  -------
//  This module is intentionally decoupled from the physics
//  core.  Its only job is to read a structured text file and
//  populate a RunParams object.  It has no knowledge of what
//  the parameters mean physically; that belongs to Calculator.
//
//  DESIGN PRINCIPLES
//  -----------------
//  1. The input file is human-readable, self-documenting and
//     INI-style:  key = value   (# comments, blank lines ok)
//  2. All parameters have hard-coded defaults (same as the
//     original Fortran hard-coded values) so that a minimal
//     input file only needs to specify what differs.
//  3. Unknown keys produce a WARNING (not a fatal error) so
//     that future parameters can be added without breaking
//     older input files.
//  4. After parsing, InputReader::validate() checks physical
//     consistency (e.g. xmin < xmax, valid enum values) and
//     returns a list of human-readable error strings.
//  5. The reader can also WRITE a template input file with all
//     keys documented, useful as a starting point.
//
//  USAGE IN YOUR OWN CODE
//  ----------------------
//    #include "InputReader.hh"
//    PionCloud::InputReader reader;
//    PionCloud::RunParams p = reader.read("myrun.in");
//    // or equivalently:
//    PionCloud::RunParams p;
//    reader.readInto("myrun.in", p);   // merges into existing p
//
//  The reader is completely stateless after construction and
//  is safe to reuse across multiple files / threads.
// ============================================================

#include "PhysicsParams.hh"
#include <string>
#include <vector>

namespace PionCloud {

// ----------------------------------------------------------
//  ParseWarning  -- non-fatal issue found during parsing
// ----------------------------------------------------------
struct ParseWarning {
    int         line;      // 1-based line number, 0 if N/A
    std::string message;
};

// ----------------------------------------------------------
//  ValidationError  -- physical inconsistency after parsing
// ----------------------------------------------------------
struct ValidationError {
    std::string field;
    std::string message;
};

// ----------------------------------------------------------
//  InputReader
// ----------------------------------------------------------
class InputReader {
public:
    InputReader() = default;

    // ---- primary interface ----

    /// Read file, return fully populated RunParams.
    /// Unset keys retain their default values.
    /// Throws std::runtime_error if the file cannot be opened.
    RunParams read(const std::string &filename);

    /// Merge file contents into an existing RunParams object.
    /// Only keys present in the file are overwritten.
    void readInto(const std::string &filename, RunParams &params);

    // ---- diagnostics (available after read/readInto) ----

    const std::vector<ParseWarning>    &warnings()  const { return warnings_; }
    const std::vector<ValidationError> &errors()    const { return errors_;   }
    bool hasErrors()   const { return !errors_.empty(); }
    bool hasWarnings() const { return !warnings_.empty(); }

    /// Print all warnings and errors to stderr.
    void printDiagnostics() const;

    // ---- utility ----

    /// Write a fully documented template input file to 'filename'.
    /// Useful for generating a starting point for new runs.
    static void writeTemplate(const std::string &filename,
                              const RunParams   &defaults = RunParams{});

    /// Validate a RunParams object independently of parsing.
    /// Returns list of errors (empty = OK).
    static std::vector<ValidationError> validate(const RunParams &p);

private:
    std::vector<ParseWarning>    warnings_;
    std::vector<ValidationError> errors_;

    // internal helpers
    void parseLine(const std::string &line, int lineNo,
                   RunParams &params);
    static std::string trim(const std::string &s);
    static std::string toLower(const std::string &s);
};

} // namespace PionCloud
