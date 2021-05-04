//

#include "backend.h"

#include "../../utils/llvm/utils.h"
#include "../../utils/llvm/codegen.h"
#include "../../utils/llvm/compile.h"
#include "../../utils/llvm/execute.h"
#include "../../utils/llvm/global_rename.h"
#include "../../utils/processor.h"
#include "../../utils/streams.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO.h>

#include <memory>

namespace NaCs::Seq::AWG {

static Device::Register<Backend> register_backend("awg");

struct Backend::ChannelInfo {
    std::string dev_name;
    uint32_t phys_chn;
    uint32_t chn_num;

    ChannelInfo(std::string &&dev_name, uint32_t phys_chn, uint32_t chn_num)
        : dev_name(std::move(dev_name)),
          phys_chn(phys_chn),
          chn_num(chn_num)
    {
    }
};

struct Backend::BasicSeq{
    
};

}
