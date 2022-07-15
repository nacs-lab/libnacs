/*************************************************************************
 *   Copyright (c) 2019 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "env.h"

#include "../nacs-utils/llvm/analysis.h"
#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/llvm/global_rename.h"
#include "../nacs-utils/llvm/utils.h"
#include "../nacs-utils/streams.h"

#include <llvm/ADT/SetVector.h>
#include <llvm/ADT/StringMap.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/Utils/Cloning.h>

#include <map>
#include <vector>

namespace NaCs::Seq {

namespace {
struct DefaultCGContext: LLVM::Codegen::CachedContext {
    DefaultCGContext(llvm::LLVMContext &llvm_ctx)
        : LLVM::Codegen::CachedContext(llvm_ctx)
    {
    }
private:
    bool use_extern_data() const override
    {
        return true;
    }
    void add_extern_data(llvm::StringRef name, const void *_data, size_t size) override
    {
        auto data = (const double*)_data;
        m_data[name] = std::vector<double>{data, data + size / sizeof(double)};
    }
    uintptr_t get_extern_data(llvm::StringRef name) override
    {
        auto it = m_data.find(name);
        if (it != m_data.end())
            return (uintptr_t)&it->second[0];
        return 0;
    }
    llvm::StringMap<std::vector<double>> m_data;
};
}

NACS_EXPORT_ Env::Env(std::unique_ptr<LLVM::Codegen::Context> cgctx)
    : m_llvm_mod(new llvm::Module("", cgctx->get_context())),
      m_cgctx(std::move(cgctx))
{
    m_cgctx->set_module(m_llvm_mod.get());
}

NACS_EXPORT_ Env::Env(llvm::LLVMContext &llvm_ctx)
    : Env(std::make_unique<DefaultCGContext>(llvm_ctx))
{
}

NACS_EXPORT() Env::~Env()
{
    while (m_vars) {
        m_vars->release();
    }
}

NACS_INTERNAL Var *Env::new_var()
{
    auto var = new Var(*this);
    link_var(var);
    m_varid_dirty = true;
    return var;
}

NACS_INTERNAL void Env::link_var(Var *var)
{
    if (m_vars)
        m_vars->m_prev = &var->m_next;
    var->m_prev = &m_vars;
    var->m_next = m_vars;
    m_vars = var;
}

NACS_EXPORT() Var *Env::new_const(IR::TagVal c)
{
    auto var = new_var();
    var->assign_const(c);
    return var;
}

static void check_args(llvm::ArrayRef<Arg> args, int nfreeargs)
{
    for (auto &arg: args) {
        if (arg.is_arg() && arg.get_arg() >= nfreeargs) {
            throw std::out_of_range("Argument index out of bound.");
        }
        else if (arg.is_var() && arg.get_var()->nfreeargs() > 0) {
            throw std::invalid_argument("Argument is not fixed.");
        }
    }
}

NACS_EXPORT() Var *Env::new_call(Var *func, llvm::ArrayRef<Arg> args, int nfreeargs)
{
    // Allow the callee to have 0 argument
    // to be consistent with what we allow in optimiztion
    if (func->nfreeargs() > 0 && func->nfreeargs() != (ssize_t)args.size())
        throw std::invalid_argument("Argument number mismatch");
    check_args(args, nfreeargs);
    // There are significantly fewer things to check since we assume the callee
    // is already verified.
    auto var = new_var();
    var->_assign_call(func, args, nfreeargs);
    return var;
}

static bool check_type(llvm::Type *ty)
{
    return uint8_t(LLVM::get_ir_type(ty)) > 0;
}

NACS_INTERNAL Var *Env::_new_call(llvm::Function *func, llvm::ArrayRef<Arg> args,
                                  int nfreeargs)
{
    if (func->arg_size() != args.size())
        throw std::invalid_argument("Argument number mismatch");
    check_args(args, nfreeargs);
    auto var = new_var();
    var->_assign_call(func, args, nfreeargs);
    return var;
}

NACS_EXPORT() Var *Env::new_call(llvm::Function *func, llvm::ArrayRef<Arg> args,
                                 int nfreeargs)
{
    auto fty = func->getFunctionType();
    if (!check_type(fty->getReturnType()))
        throw std::invalid_argument("Invalid return type.");
    if (func->isVarArg())
        throw std::invalid_argument("Vararg function not supported.");
    for (auto t: fty->params()) {
        if (!check_type(t)) {
            throw std::invalid_argument("Invalid argumentt type.");
        }
    }
    return _new_call(func, args, nfreeargs);
}

NACS_EXPORT() Var *Env::new_call(const IR::Function &func, llvm::ArrayRef<Arg> args,
                                 int nfreeargs)
{
    return _new_call(cg_context()->emit_function(func, "f"), args, nfreeargs);
}

NACS_EXPORT() Var *Env::new_extern(std::pair<IR::Type,uint64_t> ext)
{
    auto var = new_var();
    var->assign_extern(ext);
    return var;
}

NACS_EXPORT() void Env::compute_varid()
{
    int nvars = 0;
    for (auto var NACS_UNUSED: *this)
        nvars++;
    for (auto var: *this)
        var->m_varid = --nvars;
    m_varid_dirty = false;
    assert(nvars == 0);
}

NACS_EXPORT() void Env::compute_varuse()
{
    // Collect all roots.
    llvm::SmallVector<Var*, 16> roots;
    for (auto var: *this) {
        var->m_used = false;
        if (var->m_extern_ref > 0) {
            roots.push_back(var);
        }
    }

    struct VarUseVisitor : Visitor {
        bool previsit(Var *var)
        {
            if (var->m_used)
                return false;
            var->m_used = true;
            // If this is a externally and internally referenced variable,
            // it'll be scanned from the root so we can stop here.
            if (var->m_extern_ref > 0 || var->ref_none())
                return false;
            return true;
        }
    };

    VarUseVisitor visitor;
    for (auto var: roots) {
        if (var->ref_none())
            continue;
        scan_dfs(var, visitor);
    }

    m_varuse_dirty = false;
}

NACS_EXPORT() void Env::gc()
{
    compute_varuse();
    for (auto _var = m_vars; _var;) {
        auto var = _var;
        _var = _var->m_next;
        if (var->m_extern_ref > 0 || var->m_used)
            continue;
        var->release();
        m_varid_dirty = true;
    }
    // Remove unused variables should not affect varuse.
    assert(!m_varuse_dirty);
}

bool Env::optimize_local()
{
    assert(!m_varid_dirty);
    if (!m_vars)
        return false;
    bool changed = false;
    // 0: unhandled
    // 1: working on it
    // 2: done
    llvm::SmallVector<uint8_t, 64> var_states(m_vars->m_varid + 1, 0); // nvars
    llvm::SmallVector<std::pair<Var*,ssize_t>, 16> stack;
    ssize_t idx;
    Var *var;
    auto next = [&] {
        // Next argument
        if (++idx <= (ssize_t)var->args().size())
            return true;
        var_states[var->m_varid] = 2;
        // Done
        if (stack.empty())
            return false;
        std::tie(var, idx) = stack.pop_back_val();
        return true;
    };
    auto visit_field = [&] (Var *new_var) {
        auto &vs = var_states[new_var->m_varid];
        if (vs == 2)
            return next();
        if (!new_var->is_call()) {
            vs = 2;
            return next();
        }
        if (vs == 1)
            throw std::runtime_error("Dependency loop detected.");
        vs = 1;
        if (idx + 1 <= (ssize_t)var->args().size()) {
            stack.emplace_back(var, idx + 1);
        }
        else {
            var_states[var->m_varid] = 2;
        }
        var = new_var;
        idx = -2;
        return true;
    };
    auto iterate = [&] {
        if (idx == -2) {
            // Optimize callee
            assert(var->is_call());
            auto f = var->get_callee();
            if (!f.is_llvm)
                return visit_field(f.var);
            return next();
        }
        else if (idx == -1) {
            // Check callee
            changed |= var->inline_callee();
            return next();
        }
        else if (idx < (ssize_t)var->args().size()) {
            // Optimize arguments
            // This must be a call since there are at least 1 argument.
            assert(var->is_call());
            auto f = var->get_callee();
            if (f.is_llvm) {
                if (LLVM::Analysis::argument_unused(*f.llvm, idx)) {
                    // Unused, ignore.
                    changed |= !var->args()[idx].is_const();
                    var->set_arg(idx, Arg::create_const(false));
                    return next();
                }
            }
            else {
                if (f.var->argument_unused(idx)) {
                    // Unused, ignore.
                    changed |= !var->args()[idx].is_const();
                    var->set_arg(idx, Arg::create_const(false));
                    return next();
                }
            }
            auto &arg = var->args()[idx];
            if (!arg.is_var())
                return next();
            auto arg_var = arg.get_var();
            if (arg_var->is_const()) {
                var->set_arg(idx, Arg::create_const(arg_var->get_const()));
                changed = true;
                return next();
            }
            return visit_field(arg_var);
        }
        else {
            // Check arguments
            if (!var->is_call())
                return next();
            int max_arg_idx = -1;
            // There could be constant variable argument again
            // since the argument var optimization might have optimized
            // the variable to a constant
            // Normalize that first.
            for (size_t i = 0; i < var->args().size(); i++) {
                auto &arg = var->args()[i];
                if (arg.is_arg()) {
                    auto idx = arg.get_arg();
                    if (idx > max_arg_idx)
                        max_arg_idx = idx;
                    continue;
                }
                if (!arg.is_var())
                    continue;
                auto arg_var = arg.get_var();
                if (arg_var->is_const()) {
                    var->set_arg(i, Arg::create_const(arg_var->get_const()));
                    changed = true;
                }
                else if (auto copy = arg_var->get_assigned_var()) {
                    var->set_arg(i, Arg::create_var(copy));
                    changed = true;
                }
            }
            // Allow removing unused trailing arguments.
            if (max_arg_idx + 1 < var->nfreeargs()) {
                var->m_n_freeargs = max_arg_idx + 1;
                changed = true;
            }
            // Now all of the non-const arguments are used by the function
            // and all of the constant arguments are stored as constant `Arg`
            // (instead of `Arg` of constant `Var`)
            // We can execute the function iff all the arguments are constant.
            if (var->optimize_call()) {
                changed = true;
                if (var->is_call()) {
                    // If `Var::optimize_call` returns `true`,
                    // we have changed the LLVM function
                    // and should check the arguments property again
                    // if it is still a call.
                    // See the `test_return_arg_indirect` test case.
                    idx = -1;
                    return next();
                }
            }
            return next();
        }
    };
    for (auto root: *this) {
        if (root->m_extern_ref <= 0)
            continue;
        if (!root->is_call()) {
            var_states[root->m_varid] = 2;
            continue;
        }
        auto &vs = var_states[root->m_varid];
        // We are just starting a cycle so there shouldn't be nothing in progress.
        assert(vs != 1);
        // Already processed
        if (vs == 2)
            continue;
        var = root;
        idx = -2;
        vs = 1;
        while (iterate()) {
        }
    }
    return changed;
}

#ifndef NDEBUG
static void assert_topological_order(Var *var)
{
    if (!var->is_call())
        return;
    auto self_id = var->varid();
    auto check = [&] (Var *child) {
        assert(child->varid() < self_id);
    };
    auto f = var->get_callee();
    if (!f.is_llvm)
        check(f.var);
    for (auto arg: var->args()) {
        if (!arg.is_var())
            continue;
        check(arg.get_var());
    }
}

static void assert_topological_order(Env &env)
{
    for (auto var: env) {
        assert_topological_order(var);
    }
}
#endif

NACS_INTERNAL bool Env::optimize_global()
{
    assert(!m_varid_dirty);
    if (!m_vars)
        return false;
#ifndef NDEBUG
    // The variables should be in topological order by construction
    // and the optimizations should keep it that way.
    assert_topological_order(*this);
#endif
    // First we need to compute the set of root variables,
    // as well as the roots all other variables are going to be inlined into.
    // The roots are the ones that are either used by an external user or
    // more than one other roots.

    // In order to compute this, we can start with the topological order of the DAG.
    // The root for a variable is then only determined by all the variables
    // in front of it (no one after it uses it)
    // and we can assign the root for each variable by,
    //
    // 1. If the variable has external use, it is its own root and we are done.
    // 2. This function is always called after `gc()` so the variable now must have
    //    at least one internal use. If all of the users have the same root, the
    //    variable will have the same root.
    // 3. Otherwise, the variable is its own root.
    //
    // Since we don't have a way to iterate the user,
    // we'll do a slightly tweaked version by marking the usee as we go.
    unsigned nvars = m_vars->m_varid + 1;
    unsigned varid = nvars;
    llvm::SmallVector<Var*, 32> roots(nvars, nullptr);
    for (auto var: *this) {
        varid--;
        assert(var->m_varid == (int)varid);
        assume(var->m_varid == (int)varid);
        Var *root;
        assert(var->m_extern_ref > 0 || !var->is_const());
        if (var->m_extern_ref > 0) {
            root = var;
            roots[varid] = root;
        }
        else {
            root = roots[varid];
            // This is after `gc()` so the variable must have a user already.
            assert(root);
        }
        if (!var->is_call())
            continue;
        // Mark all usee's.
        auto visit_field = [&] (Var *new_var) {
            if (!roots[new_var->m_varid]) {
                // We have not found a user yet, use the current root
                roots[new_var->m_varid] = root;
                return;
            }
            // The root is the same as us or the root is the variable itself
            // no need to change.
            if (roots[new_var->m_varid] == new_var || roots[new_var->m_varid] == root)
                return;
            // The root disagrees with us and is not already the variable itself.
            // The variable is used by more than one root and we need to make it a root.
            roots[new_var->m_varid] = new_var;
        };
        auto f = var->get_callee();
        if (!f.is_llvm)
            visit_field(f.var);
        for (auto &arg: var->args()) {
            if (arg.is_var()) {
                visit_field(arg.get_var());
            }
        }
    }
    // Now we've assigned each variable to a root, it's time to inline all variables
    // to the corresponding root.
    varid = nvars;
    bool changed = false;
    llvm::SmallVector<bool, 32> visited(nvars, false);
    llvm::SmallVector<Var*, 16> branches;
    llvm::SetVector<Var*, llvm::SmallVector<Var*, 16>> inputs;
    std::map<int,llvm::Type*> args;
    for (auto var: *this) {
        varid--;
        assert(var->m_varid == (int)varid);
        assume(var->m_varid == (int)varid);
        auto root = roots[varid];
        if (root != var) {
            // `var` is ordered so a non-root variable should be visited already.
            // We don't mark all external variables as visited so allow that too.
            assert(visited[varid] || !var->is_call());
            continue;
        }
        if (auto copy = root->get_assigned_var()) {
            if (roots[copy->m_varid] != root) {
                // This is a copy of another root.
                // Since the copy should have no other uses
                // (the user of it will use the copied variable instead)
                // this must be an external use and we can just let the user
                // fix the reference later.
                continue;
            }
            visited[copy->m_varid] = true;
            changed = true;
            // The copied variable belongs to us
            // so it has no extern use and must not be constant.
            assert(!copy->is_const());
            // This is a variable that has a single use, copy it.
            if (copy->is_extern()) {
                root->assign_extern(copy->get_extern());
                continue;
            }
            assert(copy->is_call());
            assert(copy->get_callee().is_llvm);
            root->_assign_call(copy->get_callee().llvm, copy->args(),
                               copy->nfreeargs());
            // Proceed normally for the new call.
        }
        else if (!root->is_call()) {
            // No need to mutate `visited` since no one is going to read this...
            continue;
        }
        assert(root->get_callee().is_llvm);
        auto visit_field = [&] (Var *new_var) {
            if (roots[new_var->m_varid] != root || !new_var->is_call()) {
                inputs.insert(new_var);
                return;
            }
            if (visited[new_var->m_varid])
                return;
            visited[new_var->m_varid] = true;
            branches.push_back(new_var);
        };
        branches.clear();
        inputs.clear();
        args.clear();
        branches.push_back(root);
        auto oldf = root->get_callee().llvm;
        auto oldfty = oldf->getFunctionType();
        // For each roots, collect all the branches that belongs to the same root
        // as well as all the leaves that does not have the same root.
        // We'll later need to inline all the branches
        // and put all the leaves as arguments.
        for (size_t scan_idx = 0; scan_idx < branches.size(); scan_idx++) {
            auto var = branches[scan_idx];
            if (!var->is_call())
                continue;
            auto f = var->get_callee();
            if (!f.is_llvm) {
                // Everything else should have been optimized out.
                assert(var->args().size() == 0);
                visit_field(f.var);
                continue;
            }
            auto nargs = var->args().size();
            for (unsigned i = 0; i < nargs; i++) {
                auto &arg = var->args()[i];
                if (var == root && arg.is_arg()) {
                    // There should be no duplicated argument
                    auto res = args.emplace(arg.get_arg(), oldfty->getParamType(i));
                    assert(res.second);
                    (void)res;
                    continue;
                }
                // Arguments cannot have free argument
                // and constants should have been inlined
                assert(arg.is_var());
                visit_field(arg.get_var());
            }
        }
        // Nothing to inline
        if (branches.size() == 1)
            continue;
        changed = true;
        // Sort in topological order for correct (reverse) evaluation order.
        std::sort(branches.begin() + 1, branches.end(), [&] (Var *a, Var *b) {
            return a->m_varid > b->m_varid;
        });
#ifndef NDEBUG
        std::is_sorted(branches.begin(), branches.end(), [&] (Var *a, Var *b) {
            return a->m_varid > b->m_varid;
        });
        assert(branches.front() == root);
#endif
        // Now the branches are in the reverse order that we'll evaluate them
        // 1. Create the function.
        assert(!oldfty->isVarArg());
        auto rt = oldfty->getReturnType();
        auto llvm_mod = llvm_module();
        auto &llvm_ctx = llvm_mod->getContext();
        llvm::SmallVector<llvm::Type*, 8> fsig;
        decltype(root->m_args) newargs;
        for (auto arg: args) {
            fsig.push_back(arg.second);
            newargs.push_back(Arg::create_arg(arg.first));
        }
        for (auto var: inputs) {
            fsig.push_back(m_cgctx->llvm_argty(var->type()));
            newargs.push_back(Arg::create_var(var));
        }
        auto fty = llvm::FunctionType::get(rt, fsig, false);
        auto f = llvm::Function::Create(fty, llvm::GlobalValue::ExternalLinkage,
                                        oldf->getName() + ".g", llvm_mod);
        f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        f->setAttributes(llvm::AttributeList::get(llvm_ctx,
                                                  LLVM::getFnAttrs(*oldf), {}, {}));
        auto *b0 = llvm::BasicBlock::Create(llvm_ctx, "top", f);
        llvm::IRBuilder<> builder(b0);
        std::map<Var*,llvm::Value*> value_map;
        llvm::SmallVector<llvm::CallInst*, 16> calls;
        // 2. Emit each variables
        {
            unsigned argi = args.size();
            for (auto var: inputs) {
                value_map[var] = f->arg_begin() + (argi++);
            }
        }
        // This loop does not run on the first element,
        // which is `root` and we deal with it separately at the end.
        for (size_t i = branches.size() - 1; i > 0; i--) {
            auto var = branches[i];
            assert(var->is_call());
            auto f = var->get_callee();
            if (!f.is_llvm) {
                auto v = value_map[f.var];
                assert(v);
                value_map[var] = v;
                continue;
            }
            llvm::SmallVector<llvm::Value*, 8> call_args;
            for (auto &arg: var->args()) {
                assert(arg.is_var());
                auto v = value_map[arg.get_var()];
                assert(v);
                auto ty = f.llvm->getArg(call_args.size())->getType();
                // Hard code this for now (assume T_i8 are booleans)
                if (ty == m_cgctx->T_i8)
                    v = LLVM::convert_scalar(builder, m_cgctx->T_bool, v);
                call_args.push_back(LLVM::convert_scalar(builder, ty, v));
            }
            auto call = builder.CreateCall(f.llvm, call_args);
            calls.push_back(call);
            value_map[var] = call;
        }
        llvm::CallInst *res;
        {
            std::map<unsigned,unsigned> arg_map;
            llvm::SmallVector<llvm::Value*, 8> call_args;
            unsigned argi = 0;
            for (auto arg: args)
                arg_map[arg.first] = argi++;
            for (auto &arg: root->args()) {
                auto ty = oldf->getArg(call_args.size())->getType();
                llvm::Value *v;
                if (arg.is_var()) {
                    v = value_map[arg.get_var()];
                }
                else {
                    assert(arg.is_arg());
                    v = f->arg_begin() + arg_map[arg.get_arg()];
                }
                assert(v);
                // Hard code this for now (assume T_i8 are booleans)
                if (ty == m_cgctx->T_i8)
                    v = LLVM::convert_scalar(builder, m_cgctx->T_bool, v);
                call_args.push_back(LLVM::convert_scalar(builder, ty, v));
            }
            res = builder.CreateCall(oldf, call_args);
        }
        calls.push_back(res);
        builder.CreateRet(res);
        // Inline function calls
        llvm::InlineFunctionInfo IFI;
        for (auto call: calls) {
#if LLVM_VERSION_MAJOR >= 11
            llvm::InlineFunction(*call, IFI);
#else
            llvm::InlineFunction(call, IFI);
#endif
        }
        root->_assign_call(f, newargs, root->nfreeargs());
    }
    return changed;
}

NACS_EXPORT() bool Env::ensure_sorted()
{
    if (!m_sort_dirty)
        return false;

    for (auto var: *this)
        var->m_used = false;

    struct SortVisitor : Visitor {
        bool previsit(Var *var)
        {
            if (var->m_used)
                return false;
            var->m_used = true;
            if (var->ref_none()) {
                add_to_order(var);
                return false;
            }
            return true;
        }
        void postvisit(Var *var)
        {
            add_to_order(var);
        }
        llvm::SmallVector<Var*, 16> order;
        void add_to_order(Var *var)
        {
            assert(var->m_used);
            var->m_varid = order.size();
#ifndef NDEBUG
            assert_topological_order(var);
#endif
            order.push_back(var);
        }
    };

    // Mark varid as clean so that the topological order check can work.
    m_varid_dirty = false;
    SortVisitor visitor;

    for (auto var: *this) {
        if (var->m_used)
            continue;
        var->m_used = true;
        if (var->ref_none()) {
            visitor.add_to_order(var);
            continue;
        }
        scan_dfs(var, visitor);
    }
    auto &order = visitor.order;
    auto nvar = order.size();
    m_vars = nullptr;
    for (unsigned i = 0; i < nvar; i++) {
        auto var = order[i];
        link_var(var);
        assert((unsigned)var->m_varid == i);
    }
    m_varuse_dirty = true; // We've misused m_used for something else....
    m_sort_dirty = false;
    return true;
}

NACS_EXPORT() void Env::optimize()
{
    if (!m_vars)
        return;

    ensure_sorted();

    auto gc_compute_id = [&] {
        gc();
        if (m_varid_dirty) {
            compute_varid();
        }
    };

    // Initialize variable ID eagerly so the optimization functions don't need to check
    if (m_varid_dirty)
        compute_varid();
    // Local optimization and global optimization should be complete on their own
    // and won't expose more optimization opportunities
    // for themselves when called back to back.
    // However, they could expose optimization opportunities for each other
    // so we should call the other one if the one of them made a change.
    // Additionally, we should call each one (along with `gc()`) at least once
    // to capture all the initial optimization chances.
    optimize_local();
    // `optimize_global` assumes all nodes are used so we need to run `gc()` once first.
    gc_compute_id();
    while (true) {
        if (!optimize_global())
            break;
        gc_compute_id();
        if (!optimize_local())
            break;
        gc_compute_id();
    }
    finalize_vars();
}

void Env::finalize_vars()
{
    // We don't use any LLVM global info for optimization so we only need to do the LLVM
    // global DCE once.
    for (auto &go: llvm_module()->global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    auto type_matches = [&] (Var *var) {
        auto f = var->get_callee();
        assert(f.is_llvm);
        auto args = var->args();
        uint32_t nargs = args.size();
        auto fty = f.llvm->getFunctionType();
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            if (arg.is_arg())
                continue;
            // Constants should have been optimized out.
            assert(arg.is_var());
            if (fty->getParamType(i) != arg.get_var()->llvm_type()) {
                return false;
            }
        }
        return true;
    };
    // Make sure the argument type matches the LLVM type
    for (auto var: *this) {
        if (!var->is_call())
            continue;
        auto f = var->get_callee();
        // There can be assignments
        if (!f.is_llvm) {
            assert(var->args().empty());
            continue;
        }
        else if (type_matches(var)) {
            f.llvm->setLinkage(llvm::GlobalValue::ExternalLinkage);
            f.llvm->setVisibility(llvm::GlobalValue::ProtectedVisibility);
            continue;
        }
        auto args = var->args();
        uint32_t nargs = args.size();
        llvm::SmallVector<llvm::Type*,8> arg_types(nargs);
        auto oldfty = f.llvm->getFunctionType();
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            if (arg.is_arg()) {
                arg_types[i] = oldfty->getParamType(i);
            }
            else {
                // Constants should have been optimized out.
                assert(arg.is_var());
                arg_types[i] = arg.get_var()->llvm_type();
            }
        }
        auto rt = oldfty->getReturnType();
        auto llvm_mod = llvm_module();
        auto &llvm_ctx = llvm_mod->getContext();
        auto fty = llvm::FunctionType::get(rt, arg_types, false);
        auto newf = llvm::Function::Create(fty, llvm::GlobalValue::ExternalLinkage,
                                           f.llvm->getName() + ".t", llvm_mod);
        newf->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        newf->setAttributes(llvm::AttributeList::get(
                                llvm_ctx, LLVM::getFnAttrs(*f.llvm), {}, {}));
        auto *b0 = llvm::BasicBlock::Create(llvm_ctx, "top", newf);
        llvm::IRBuilder<> builder(b0);
        llvm::SmallVector<llvm::Value*,8> call_args(nargs);
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            llvm::Value *larg = newf->getArg(i);
            if (!arg.is_arg()) {
                // Constants should have been optimized out.
                assert(arg.is_var());
                if (oldfty->getParamType(i) != larg->getType()) {
                    if (arg.get_var()->type() == IR::Type::Bool ||
                        oldfty->getParamType(i) == m_cgctx->T_i8)
                        larg = LLVM::convert_scalar(builder, m_cgctx->T_bool, larg);
                    larg = LLVM::convert_scalar(builder, oldfty->getParamType(i), larg);
                }
            }
            call_args[i] = larg;
        }
        auto call = builder.CreateCall(f.llvm, call_args);
        builder.CreateRet(call);
        llvm::InlineFunctionInfo IFI;
#if LLVM_VERSION_MAJOR >= 11
        llvm::InlineFunction(*call, IFI);
#else
        llvm::InlineFunction(call, IFI);
#endif
        var->_assign_call(newf, args, var->nfreeargs());
    }
    llvm::legacy::PassManager PM;
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.add(llvm::createGlobalDCEPass());
    // Shorten all global names since we don't care what they are
    // and this should slightly reduce the compiled binary size.
    PM.add(LLVM::createGlobalRenamePass());
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif

    PM.run(*llvm_module());
}

NACS_EXPORT() void Env::print(std::ostream &stm, bool sortvar) const
{
    if (sortvar)
        const_cast<Env*>(this)->ensure_sorted();
    llvm::SmallVector<Var*, 32> vars(begin(), end());
    if (vars.empty()) {
        stm << "Variables: <empty>" << std::endl;
    }
    else {
        stm << "Variables: <" << vars.size() << ">" << std::endl;
    }
    for (auto it = vars.rbegin(), end = vars.rend(); it != end; ++it) {
        stm << "  ";
        (*it)->print(stm);
        stm << std::endl;
    }
    stm << std::endl;
    if (!llvm_module()) {
        stm << "LLVM: <null>" << std::endl;
    }
    else {
        stm << "LLVM:";
        string_ostream sstm;
        llvm::raw_os_ostream lstm(sstm);
        llvm_module()->print(lstm, nullptr);
        auto str = sstm.get_buf();
        auto p = &str[0];
        auto end = p + str.size();
        for (; p < end && *p == '\n'; p++) {
        }
        if (p < end) {
            stm << std::endl;
            while (p < end) {
                auto lend = (char*)memchr(p, '\n', end - p);
                stm.write("  ", 2);
                if (!lend) {
                    stm.write(p, end - p);
                    stm.put('\n');
                    break;
                }
                stm.write(p, lend + 1 - p);
                p = lend + 1;
            }
        }
        else {
            stm << " <empty>" << std::endl;
        }
    }
}

}
