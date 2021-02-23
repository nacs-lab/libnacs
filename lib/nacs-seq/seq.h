/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_SEQ_SEQ_H__
#define __NACS_SEQ_SEQ_H__

#include "basic_seq.h"

#include <list>
#include <unordered_map>

#include <llvm/ADT/StringMap.h>

namespace NaCs::Seq {

class Seq {
public:
    Seq(std::unique_ptr<LLVM::Codegen::Context> cgctx);
    Seq(llvm::LLVMContext &llvm_ctx);
    ~Seq();
    Env &env()
    {
        return m_env;
    }
    const Env &env() const
    {
        return m_env;
    }
    BasicSeq *add_basicseq(uint32_t id);
    uint32_t get_chn_id(llvm::StringRef name, bool create=false);
    const std::string &get_chn_name(uint32_t id) const;
    llvm::ArrayRef<std::string> get_chn_names() const
    {
        return m_chnnames;
    }
    const std::list<BasicSeq> &get_basicseqs() const
    {
        return m_seqs;
    }
    const std::map<uint32_t,double> &get_defvals() const
    {
        return m_defval;
    }

    void set_defval(uint32_t chn, double val);
    double defval(uint32_t chn);

    Var *get_const(IR::TagVal c);
    Var *get_call(const IR::Function &func, llvm::ArrayRef<Arg> args, int nfreeargs);
    // `id` is expected to be a counter from 0.
    // This will be used as the index into an array of sequence variables.
    // If different `typ`s are used with the same `id`,
    // only the `typ` for the first call will be used.
    Var *get_slot(IR::Type typ, uint32_t id);
    llvm::ArrayRef<Var::Ref> get_slots() const
    {
        return m_slot_vars;
    }

    void optimize();
    void prepare();
    void print(std::ostream &stm) const;

private:
    bool optimize_cfg();
    bool optimize_chn(uint32_t chn);

    // This should be the first member so that this is destroyed last,
    // before everything else that refers to this.
    Env m_env;
    std::list<BasicSeq> m_seqs;
    std::unordered_map<double,Var*> m_float_cache;
    std::map<uint32_t,double> m_defval;
    std::vector<Var::Ref> m_slot_vars;
    llvm::StringMap<uint32_t> m_chnids;
    std::vector<std::string> m_chnnames;
};

static inline std::ostream &operator<<(std::ostream &stm, const Seq &seq)
{
    seq.print(stm);
    return stm;
}

}

#endif
