#pragma once
// ============================================================
//  Models.hh   (new in v1.4)
//  Pluggable physics ingredients of the pion-cloud calculation.
//
//  WHY
//  The tagged structure function factorises into three independent
//  ingredients:
//      F2^(MB)(x; cuts) = int [ f_MB(y,kT) ] x [ F2^M(x/y, Q2) ] over the
//                         region allowed by the cuts (kinematics).
//  Up to v1.3 the first two were hard-wired with switch statements in
//  several places, so a new model (chiral-EFT splitting function, kaon
//  cloud, a modern pion PDF, ...) meant editing the integration code.
//  Now the Calculator only talks to these two interfaces:
//
//    SplittingModel : f(y, kT)        per d(kT^2), see SplittingFunctions.hh
//    PionSFModel    : F2(x_pi, Q2)    meson structure function
//
//  The built-in models below wrap the original functions and give
//  bit-identical results (verified on all validation runs).
//
//  HOW A THEORIST ADDS A MODEL (no change to PionCloud's sources):
//
//    #include "Models.hh"
//    struct MyChiralPiN : PionCloud::SplittingModel {
//        double f(double y, double kT) const override { ... }
//        std::string name() const override { return "piN_chiralEFT"; }
//    };
//    // either inject it directly:
//    PionCloud::Calculator calc(params, std::make_shared<MyChiralPiN>());
//    // or register it by name, then select it from an input file with
//    //   splitting_model = piN_chiralEFT
//    PionCloud::ModelRegistry::addSplitting("piN_chiralEFT",
//        [](const PionCloud::RunParams &p){ return std::make_shared<MyChiralPiN>(); });
//
//  The same works for PionSFModel (e.g. a table or LHAPDF set):
//  ModelRegistry::addPionSF("name", factory) / input key pion_sf_model.
//
//  Thread safety: models must be const-callable from several threads
//  if the caller parallelises (the built-ins are pure functions).
// ============================================================

#include "PhysicsParams.hh"
#include "SplittingFunctions.hh"
#include "GRV.hh"

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PionCloud {

// ----------------------------------------------------------
//  Interfaces
// ----------------------------------------------------------
struct SplittingModel {
    virtual ~SplittingModel() = default;
    // Light-cone probability density per d(kT^2): y = meson fraction, kT [GeV].
    virtual double f(double y, double kT) const = 0;
    virtual std::string name() const = 0;
};

struct PionSFModel {
    virtual ~PionSFModel() = default;
    // Meson structure function F2^M(x_M, Q2).
    virtual double F2(double xpi, double Q2) const = 0;
    virtual std::string name() const = 0;
};

// ----------------------------------------------------------
//  Built-in splitting functions (the Fortran ports)
// ----------------------------------------------------------
class BuiltinSplitting : public SplittingModel {
public:
    explicit BuiltinSplitting(const RunParams &p)
        : flag_(p.flag), L_(p.L), Ld_(p.Ld), typ_(p.typ), dis_(p.dis) {}
    double f(double y, double kT) const override {
        switch (flag_) {
        case FLAG_PION:      return fypiN   (y, kT, L_,  typ_, dis_);
        case FLAG_RHO:       return f_rhoN  (y, kT, L_,  typ_, dis_);
        case FLAG_PI_DELTA:  return fypiD   (y, kT, Ld_, typ_, dis_);
        case FLAG_RHO_DELTA: return f_RhoDel(y, kT, Ld_, typ_, dis_);
        }
        throw std::invalid_argument("BuiltinSplitting: unknown ContributionFlag");
    }
    std::string name() const override {
        switch (flag_) {
        case FLAG_PION:      return "pion";
        case FLAG_RHO:       return "rho";
        case FLAG_PI_DELTA:  return "pi_delta";
        case FLAG_RHO_DELTA: return "rho_delta";
        }
        return "unknown";
    }
private:
    ContributionFlag flag_;
    double L_, Ld_;
    FormFactorType typ_;
    DisChannel dis_;
};

// ----------------------------------------------------------
//  Built-in pion structure functions (GRV92 / GRS99, LO)
// ----------------------------------------------------------
class BuiltinPionSF : public PionSFModel {
public:
    explicit BuiltinPionSF(const RunParams &p)
        : grv_(p.grv_version), strange_(p.pion_strange) {}
    double F2(double xpi, double Q2) const override {
        return (grv_ == GRV_99) ? F2piGRV99(xpi, Q2, strange_)
                                : F2piGRV  (xpi, Q2, strange_);
    }
    std::string name() const override {
        return std::string(grv_ == GRV_99 ? "GRS99" : "GRV92") +
               (strange_ ? "+s" : "");
    }
private:
    GRVVersion grv_;
    bool strange_;
};

// ----------------------------------------------------------
//  Registry: models selectable by name from input files
// ----------------------------------------------------------
class ModelRegistry {
public:
    using SplitFactory = std::function<std::shared_ptr<const SplittingModel>(const RunParams &)>;
    using SFFactory    = std::function<std::shared_ptr<const PionSFModel>(const RunParams &)>;

    static void addSplitting(const std::string &name, SplitFactory f) { splitMap()[name] = std::move(f); }
    static void addPionSF   (const std::string &name, SFFactory    f) { sfMap()[name]    = std::move(f); }

    // Empty name -> built-in model chosen by p.flag / p.grv_version.
    static std::shared_ptr<const SplittingModel> makeSplitting(const RunParams &p) {
        if (p.splitting_model.empty()) return std::make_shared<BuiltinSplitting>(p);
        auto it = splitMap().find(p.splitting_model);
        if (it == splitMap().end())
            throw std::invalid_argument("ModelRegistry: no splitting model named '" +
                                        p.splitting_model + "' is registered");
        return it->second(p);
    }
    static std::shared_ptr<const PionSFModel> makePionSF(const RunParams &p) {
        if (p.pion_sf_model.empty()) return std::make_shared<BuiltinPionSF>(p);
        auto it = sfMap().find(p.pion_sf_model);
        if (it == sfMap().end())
            throw std::invalid_argument("ModelRegistry: no pion SF model named '" +
                                        p.pion_sf_model + "' is registered");
        return it->second(p);
    }
    static std::vector<std::string> splittingNames() { return keys(splitMap()); }
    static std::vector<std::string> pionSFNames()    { return keys(sfMap()); }

private:
    static std::map<std::string, SplitFactory> &splitMap() { static std::map<std::string, SplitFactory> m; return m; }
    static std::map<std::string, SFFactory>    &sfMap()    { static std::map<std::string, SFFactory> m; return m; }
    template <class M> static std::vector<std::string> keys(const M &m) {
        std::vector<std::string> v; for (auto &kv : m) v.push_back(kv.first); return v;
    }
};

} // namespace PionCloud
