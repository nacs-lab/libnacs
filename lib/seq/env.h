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

#ifndef __NACS_SEQ_ENV_H__
#define __NACS_SEQ_ENV_H__

#include <nacs-utils/utils.h>
#include <nacs-utils/ir.h>

#include <llvm/IR/Module.h>

#include <iterator>
#include <memory>

namespace NaCs::LLVM::Codegen {
class Context;
}

namespace NaCs::Seq {

class Env;
class Var;

class Arg {
    enum ArgType : uint8_t {
        Const,
        Var_,
        Arg_,
    };
public:
    static Arg create_const(IR::TagVal c)
    {
        Arg a;
        a.m_valtyp = c.typ;
        a.m_argtyp = Const;
        a.m_const = c.val;
        return a;
    }
    static Arg create_var(Var *var)
    {
        Arg a;
        a.m_argtyp = Var_;
        a.m_var = var;
        return a;
    }
    static Arg create_arg(int arg)
    {
        Arg a;
        a.m_argtyp = Arg_;
        a.m_arg = arg;
        return a;
    }
    bool is_const() const
    {
        return m_argtyp == Const;
    }
    IR::TagVal get_const() const
    {
        return IR::TagVal(m_valtyp, m_const);
    }
    bool is_var() const
    {
        return m_argtyp == Var_;
    }
    Var *get_var() const
    {
        return assume(m_var);
    }
    bool is_arg() const
    {
        return m_argtyp == Arg_;
    }
    int get_arg() const
    {
        return m_arg;
    }
    // C++20 operator<=>
    int compare(const Arg &other) const
    {
        if (m_argtyp != other.m_argtyp)
            return int(m_argtyp) - int(other.m_argtyp);
        switch (m_argtyp) {
        case ArgType::Const:
            if (m_valtyp != other.m_valtyp)
                return int(m_valtyp) - int(other.m_valtyp);
            switch (m_valtyp) {
            case IR::Type::Bool:
                return int(m_const.b) - int(other.m_const.b);
            case IR::Type::Int32:
                return m_const.i32 - other.m_const.i32;
            case IR::Type::Float64:
                if (m_const.f64 == other.m_const.f64)
                    return 0;
                return m_const.f64 < other.m_const.f64 ? -1 : 1;
            default:
                return 0;
            }
        case ArgType::Var_:
            return intptr_t(m_var) - intptr_t(other.m_var);
        case ArgType::Arg_:
            return m_arg - other.m_arg;
        default:
            return 0;
        }
    }

private:
    // Keep the value type out of the union so that the two single byte members
    // share the same alignment and the size of the whole struct can be 16.
    IR::Type m_valtyp;
    ArgType m_argtyp;
    union {
        IR::GenVal m_const;
        Var *m_var;
        int m_arg;
    };
};
static_assert(sizeof(Arg) == 8 * 2, "Check alignment.");

std::ostream &operator<<(std::ostream &stm, const Arg &arg);

class Var {
    // The user cannot mutate this freely.
    // Each variable can only refer to a variable added before it.
    // (i.e. the linked list in `Env` natually forms a topological order of the variables)
    struct Unrefer {
        void operator()(Var *var) const
        {
            var->unref();
        }
    };

public:
    struct Func {
        bool is_llvm;
        union {
            Var *var;
            llvm::Function *llvm;
        };
    };
    // Use unique pointer since I'm too lazy to write my own.
    // The only difference is that copy is not allowed
    // but we could use `->ref()` most of the time.
    using Ref = std::unique_ptr<Var,Unrefer>;

    Ref ref()
    {
        if (!m_extern_ref)
            varuse_dirty();
        m_extern_ref++;
        return Ref(this);
    }

    Env &env() const
    {
        return m_env;
    }
    llvm::ArrayRef<Arg> args() const
    {
        return m_args;
    }
    bool valid() const
    {
        return m_in_env;
    }
    uint32_t extern_used() const
    {
        return m_extern_ref;
    }
    bool is_const() const
    {
        return m_is_const;
    }
    IR::TagVal get_const() const
    {
        return m_value.constant;
    }
    bool is_extern() const
    {
        return m_is_extern;
    }
    std::pair<IR::Type,uint64_t> get_extern() const
    {
        return m_value.ext;
    }
    int nfreeargs() const
    {
        return m_n_freeargs;
    }
    bool is_call() const
    {
        return !is_const() && !is_extern();
    }
    Func get_callee() const
    {
        return m_value.func;
    }
    IR::Type type() const;
    int varid() const;
    bool used(bool ext) const;
    void print(std::ostream &stm, bool newline=false) const;
    bool argument_unused(int idx) const;
    // Return the variable that this one is a copy of
    // (i.e. zero argument call of another variable)
    Var *get_assigned_var() const
    {
        if (is_call() && !get_callee().is_llvm && nfreeargs() == 0)
            return get_callee().var;
        return nullptr;
    }

private:
    // The following functions are only used in optimizations or to create new variables
    // and does not need to explicitly invalidate varuse.
    void set_arg(int i, Arg arg)
    {
        m_args[i] = arg;
    }
    void assign_const(IR::TagVal c)
    {
        m_is_const = true;
        m_is_extern = false;
        m_n_freeargs = 0;
        m_value.constant = c;
        m_args.clear();
    }
    void assign_call(Var *func, llvm::ArrayRef<Arg> args, int nfreeargs=0);
    void assign_call(llvm::Function *func, llvm::ArrayRef<Arg> args, int nfreeargs=0);
    // Special case for assigning a copy of a variable.
    void assign_var(Var *var)
    {
        m_args.clear();
        m_n_freeargs = 0;
        m_is_const = false;
        m_is_extern = false;
        m_value.func.is_llvm = false;
        m_value.func.var = var;
    }
    void assign_extern(std::pair<IR::Type,uint64_t> ext)
    {
        m_is_const = false;
        m_is_extern = true;
        m_n_freeargs = 0;
        m_value.ext = ext;
        m_args.clear();
    }
    void fill_args(llvm::ArrayRef<Arg> args, int nfreeargs);

    bool inline_callee();
    bool optimize_call();
    void optimize_llvmf(llvm::Function *f);
    // Remove unnecessary arguments:
    // * Inline constant arguments (this should also remove unused ones)
    // * Remove duplicated arguments
    bool reduce_args();

    Var(Env &env)
        : m_env{env}
    {
    }
    ~Var()
    {
        if (m_in_env) {
            unlink();
        }
    }
    void unlink()
    {
        if (m_next)
            m_next->m_prev = m_prev;
        *m_prev = m_next;
    }
    void unref()
    {
        // Be **VERY** careful since this function can `delete this`.
        if (--m_extern_ref > 0)
            return;
        if (m_in_env) {
            varuse_dirty();
            return;
        }
        delete this;
    }
    void release()
    {
        if (m_extern_ref) {
            m_in_env = false;
            unlink();
        }
        else {
            delete this;
        }
    }
    void varuse_dirty();

    Env &m_env;
    Var **m_prev{nullptr};
    Var *m_next{nullptr};
    // If `m_extern_ref > 0`, the var is used externally
    // and shouldn't be optimized out completely.
    // (Also it won't be freed when env is freed).
    uint32_t m_extern_ref{0};
    // Variable ID. Lazily computed from the list of variable in the env.
    int m_varid;
    // Whether this variable is used by other variables that are used by external reference.
    // This will be `false` if `m_extern_ref > 0` but no other variable uses this.
    bool m_used;

    // If `m_in_env` `false` the var is not associate with a env anymore.
    bool m_in_env{true};

    bool m_is_const{false};
    bool m_is_extern{false};
    int m_n_freeargs{0};
    union {
        IR::TagVal constant;
        Func func;
        std::pair<IR::Type,uint64_t> ext{};
    } m_value;
    llvm::SmallVector<Arg, 4> m_args{};
    friend class Env;
    friend std::ostream &operator<<(std::ostream &stm, const Arg &arg);
};

static inline std::ostream &operator<<(std::ostream &stm, const Var &var)
{
    var.print(stm);
    return stm;
}

class Env {
public:
    class iterator { // input iterator
    public:
        iterator(Var *var)
            : m_var(var)
        {
        }
        bool operator==(const iterator &other) const
        {
            return m_var == other.m_var;
        }
        bool operator!=(const iterator &other) const
        {
            return m_var != other.m_var;
        }
        iterator &operator++()
        {
            m_var = m_var->m_next;
            return *this;
        }
        iterator operator++(int)
        {
            auto copy = *this;
            this->operator++();
            return copy;
        }
        Var *operator*() const
        {
            return m_var;
        }
        Var *const *operator->() const
        {
            return &m_var;
        }

    private:
        Var *m_var;
    };

    Env(std::unique_ptr<LLVM::Codegen::Context> cgctx);
    Env(llvm::LLVMContext &llvm_ctx);
    ~Env();

    iterator begin() const
    {
        return m_vars;
    }
    iterator end() const
    {
        return nullptr;
    }

    Var *new_const(IR::TagVal c);
    Var *new_call(Var *func, llvm::ArrayRef<Arg> args, int nfreeargs=0);
    Var *new_call(llvm::Function *func, llvm::ArrayRef<Arg> args, int nfreeargs=0);
    Var *new_call(const IR::Function &func, llvm::ArrayRef<Arg> args, int nfreeargs=0);
    Var *new_extern(std::pair<IR::Type,uint64_t> ext);
    llvm::Module *llvm_module() const
    {
        return m_llvm_mod.get();
    }
    LLVM::Codegen::Context *cg_context() const
    {
        return m_cgctx.get();
    }

    int num_vars() const;
    void gc();
    void optimize();
    void print(std::ostream &stm) const;

private:
    Var *_new_call(llvm::Function *func, llvm::ArrayRef<Arg> args, int nfreeargs);
    Var *new_var();
    void compute_varid();
    void compute_varuse();
    // Optimize each variables (functions) individually based only on the
    // arguments and callee without information about the global call graph.
    bool optimize_local();

    std::unique_ptr<llvm::Module> m_llvm_mod;
    std::unique_ptr<LLVM::Codegen::Context> m_cgctx;
    Var *m_vars{nullptr};

    bool m_varid_dirty{false};
    bool m_varuse_dirty{false};
    friend class Var;
};

inline int Var::varid() const
{
    if (!m_in_env)
        return -1;
    if (m_env.m_varid_dirty)
        m_env.compute_varid();
    return m_varid;
}

inline void Var::varuse_dirty()
{
    if (m_in_env) {
        m_env.m_varuse_dirty = true;
    }
}

inline bool Var::used(bool ext) const
{
    if (!m_in_env)
        return false;
    if (ext && m_extern_ref > 0)
        return true;
    if (m_env.m_varuse_dirty)
        m_env.compute_varuse();
    return m_used;
}

inline int Env::num_vars() const
{
    if (!m_vars)
        return 0;
    return m_vars->varid() + 1;
}

static inline std::ostream &operator<<(std::ostream &stm, const Env &env)
{
    env.print(stm);
    return stm;
}

}

namespace std {
template<>
class iterator_traits<NaCs::Seq::Env::iterator> {
public:
    using value_type = NaCs::Seq::Var*;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::input_iterator_tag;
};
}

#endif
