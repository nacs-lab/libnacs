/*************************************************************************
 *   Copyright (c) 2021 - 2022 Yichao Yu <yyc1992@gmail.com>             *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *                                                                       *
 *   This library is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   Lesser General Public License for more details.                     *
 *                                                                       *
 *   You should have received a copy of the GNU Lesser General Public    *
 *   License along with this library. If not,                            *
 *   see <http://www.gnu.org/licenses/>.                                 *
 *************************************************************************/

#include "event_time.h"
#include "error.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/number.h"

#include <algorithm>

namespace NaCs::Seq {

NACS_EXPORT_ EventTime::EventTime(const EventTime &other)
    : tconst(other.tconst),
      terms(other.terms.size())
{
    // Do not copy reference count
    auto nterms = other.terms.size();
    for (size_t i = 0; i < nterms; i++) {
        auto &term = terms[i];
        auto &oterm = other.terms[i];
        term.sign = oterm.sign;
        term.var.reset(oterm.var.get());
        term.id = oterm.id;
    }
}

NACS_EXPORT() bool EventTime::normalize()
{
    bool changed = false;
    for (auto &term: terms) {
        if (auto v = term.var->get_assigned_var()) {
            term.var = v->ref();
            changed = true;
        }
    }
    // Sort the terms so that we can compare them in order in `isless_terms`.
    // Don't count change in sort since that won't expose new optimization opportunity.
    std::sort(terms.begin(), terms.end(), [] (const auto &t1, const auto &t2) {
        auto c1 = t1.var->is_const();
        auto c2 = t2.var->is_const();
        if (c1 && !c2)
            return true;
        if (c2 && !c1)
            return false;
        return t1.var->varid() < t2.var->varid();
    });
    auto nterms = terms.size();
    unsigned i = 0;
    for (; i < nterms; i++) {
        auto &term = terms[i];
        if (!term.var->is_const())
            break;
        auto t = term.var->get_const().get<double>();
        if (term.sign == Sign::Pos && !(round<int64_t>(t) > 0))
            throw Error(Error::Type::EventTime, Error::EventTime::NonPosTime,
                        Error::Type::EventTime, term.id, "Positive time expected.");
        if (term.sign == Sign::NonNeg && !(t >= 0))
            throw Error(Error::Type::EventTime, Error::EventTime::NegTime,
                        Error::Type::EventTime, term.id, "Non-negative time expected.");
        tconst += round<int64_t>(t);
    }
    if (i > 0) {
        terms.erase(terms.begin(), terms.begin() + i);
        changed = true;
    }
    return changed;
}

NACS_EXPORT() Sign EventTime::isless_terms(const EventTime &vc2) const
{
    const auto &vc1 = *this;
    auto off1 = 0u;
    auto off2 = 0u;
    auto nt1 = vc1.terms.size();
    auto nt2 = vc2.terms.size();

    // If vc1 has more terms left, we definately cannot determine that it is smaller than vc2.
    auto failed = [&] { return nt1 - off1 > nt2 - off2; };
    if (failed())
        return Sign::Unknown;
    bool has_pos = false;
    while (off1 < nt1) {
        // vc2 has at least as many terms left as vc1 so it is safe to index vc2 here.
        auto &t1 = vc1.terms[off1];
        auto &t2 = vc2.terms[off2];
        if (t1.var.get() == t2.var.get()) {
            off1++;
            off2++;
            continue;
        }
        // vc2 contains a term that we may not know about in `vc1`,
        if (t2.sign == Sign::Unknown)
            return Sign::Unknown;
        if (t2.sign == Sign::Pos)
            has_pos = true; // The difference so far is positive.
        // if it is non-negetive, we can move to the next one.
        off2++;
        if (failed()) {
            return Sign::Unknown;
        }
    }
    assert(off1 == nt1);
    // Scan the extra terms in vc2. Especially since there might be unknown signs in it.
    while (off2 < nt2) {
        auto &t2 = vc2.terms[off2];
        // vc2 contains a term that we may not know about in `vc1`,
        if (t2.sign == Sign::Unknown)
            return Sign::Unknown;
        if (t2.sign == Sign::Pos)
            has_pos = true; // The difference so far is positive.
        // if it is non-negetive, we can move to the next one.
        off2++;
    }
    return has_pos ? Sign::Pos : Sign::NonNeg;
}

NACS_EXPORT() Sign EventTime::isless(const EventTime &vc2) const
{
    if (tconst > vc2.tconst)
        return Sign::Unknown;
    auto res = isless_terms(vc2);
    if (res == Sign::NonNeg && tconst < vc2.tconst)
        return Sign::Pos;
    return res;
}

NACS_EXPORT() EventTime EventTime::operator-(const EventTime &other) const
{
    assert(other.tconst <= tconst);
    assert(other.isless_terms(*this) != Sign::Unknown);
    int64_t diff_tconst = tconst - other.tconst;
    decltype(terms) diff_terms;

    const auto &vc1 = other;
    const auto &vc2 = *this;
    auto off1 = 0u;
    auto off2 = 0u;
    auto nt1 = vc1.terms.size();
    auto nt2 = vc2.terms.size();

    auto assertsize = [&] { assert(nt1 - off1 <= nt2 - off2); };
    assertsize();
    while (off1 < nt1) {
        // vc2 has at least as many terms left as vc1 so it is safe to index vc2 here.
        auto &t1 = vc1.terms[off1];
        auto &t2 = vc2.terms[off2];
        if (t1.var.get() == t2.var.get()) {
            off1++;
            off2++;
            continue;
        }
        // vc2 contains a term that we may not know about in `vc1`,
        assert(t2.sign != Sign::Unknown);
        // if it is non-negetive, we can move to the next one.
        off2++;
        assertsize();
        diff_terms.push_back({t2.sign, t2.var->ref()});
    }
    while (off2 < nt2) {
        auto &t2 = vc2.terms[off2];
        off2++;
        diff_terms.push_back({t2.sign, t2.var->ref()});
    }
    return {diff_tconst, std::move(diff_terms)};
}

NACS_EXPORT() Var *EventTime::to_var(Env &env) const
{
    unsigned narg = terms.size();
    if (narg == 0)
        return env.new_const(IR::TagVal(double(tconst)));
    // Even if we have only one term, we still need to round to integer...
    llvm::SmallVector<Arg, 8> args;
    if (tconst != 0) {
        args.push_back(Arg::create_const(IR::TagVal(double(tconst))));
        narg += 1;
    }
    auto mod = env.llvm_module();
    auto &ctx = mod->getContext();
    auto cgctx = env.cg_context();
    llvm::Function *f;
    {
        llvm::FunctionType *ftype;
        llvm::SmallVector<llvm::Type*, 8> fsig(narg, cgctx->T_f64);
        ftype = llvm::FunctionType::get(cgctx->T_f64, fsig, false);
        f = llvm::Function::Create(ftype, llvm::GlobalValue::ExternalLinkage, "+", mod);
    }
    f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
    f->addFnAttr(llvm::Attribute::AlwaysInline);
    f->addFnAttr(llvm::Attribute::Speculatable);
    f->addFnAttr(llvm::Attribute::NoRecurse);
    f->addFnAttr(llvm::Attribute::NoUnwind);
    f->addFnAttr(llvm::Attribute::ReadOnly);

    {
        auto b0 = llvm::BasicBlock::Create(ctx, "top", f);
        llvm::IRBuilder<> builder(b0);
        llvm::FastMathFlags fmf;
        // Do not set noinf since we may use divide by 0
        fmf.setNoNaNs();
        fmf.setNoSignedZeros();
        fmf.setAllowReciprocal();
#if LLVM_VERSION_MAJOR >= 7
        fmf.setAllowContract();
#endif
        builder.setFastMathFlags(fmf);
        llvm::Value *v = nullptr;
        auto T_i64 = llvm::Type::getInt64Ty(ctx);
        // lrint can be lowered directly into SSE instructions
        auto intrin = llvm::Intrinsic::getDeclaration(mod, llvm::Intrinsic::lrint,
                                                      {T_i64, cgctx->T_f64});
        for (auto &arg: f->args()) {
            llvm::Value *argv = builder.CreateCall(intrin, {&arg});
            if (!v) {
                v = argv;
                continue;
            }
            v = builder.CreateAdd(v, argv);
        }
        builder.CreateRet(builder.CreateSIToFP(v, cgctx->T_f64));
    }
    for (auto &t: terms)
        args.push_back(Arg::create_var(t.var.get()));
    return env.new_call(f, args);
}

NACS_EXPORT() void EventTime::print(std::ostream &stm, bool newline) const
{
    stm << tconst;
    for (auto &t: terms) {
        stm << " + (";
        t.var->print(stm, false, true);
        stm << ")/";
        if (t.sign == Sign::Pos) {
            stm << "p";
        } else if (t.sign == Sign::NonNeg) {
            stm << "nn";
        } else {
            stm << "u";
        }
    }
    if (newline) {
        stm << std::endl;
    }
}

NACS_EXPORT() int64_t EventTime::min_const() const
{
    for (auto &term: terms) {
        if (term.sign == Sign::Unknown) {
            return 0;
        }
    }
    return max(tconst, 0);
}

}
