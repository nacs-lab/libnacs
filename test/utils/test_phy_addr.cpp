#include "../../lib/nacs-utils/fd_utils.h"
#include "../../lib/nacs-utils/mem.h"
#include "../../lib/nacs-utils/timer.h"

#include <unistd.h>
#include <sys/mman.h>
#include <iostream>
#include <iomanip>

using namespace NaCs;

void
test_write(volatile void *ptr, uint32_t nrun)
{
    Timer timer;
    for (uint32_t i = 0;i < nrun;i++)
        Mem::write<uint64_t>(ptr, i);
    auto time = timer.elapsed();
    std::cout << "Time per write: " << std::setprecision(4)
              << double(time) / double(nrun) / 1e3
              << " us" << std::endl;
}

int
main()
{
    void *virt_addr = mmap(nullptr, page_size, PROT_READ | PROT_WRITE,
                           MAP_SHARED | MAP_ANONYMOUS | MAP_LOCKED |
                           MAP_POPULATE, -1, 0);
    void *phy_addr = getPhyAddr(virt_addr);

    void *virt_addr2 = mapPhyAddr(phy_addr, page_size);
    strcpy((char*)virt_addr2, "random string");
    printf("%s\n", (char*)virt_addr);
    strcpy((char*)virt_addr2, "random string22222");
    printf("%s\n", (char*)virt_addr);

    std::cout << "page_size: 0x" << std::hex << page_size << std::endl;
    std::cout << "phy_addr : 0x" << std::hex << phy_addr << std::endl;
    std::cout << "virt_addr: 0x" << std::hex << virt_addr << std::endl;
    std::cout << "virt_addr2: 0x" << std::hex << virt_addr2 << std::endl;
    std::cout << "phy_addr2: 0x"
              << std::hex << getPhyAddr(virt_addr2) << std::endl;

    *(volatile uint64_t*)virt_addr = 0;
    msync(virt_addr, page_size, MS_SYNC);
    msync(virt_addr2, page_size, MS_INVALIDATE);
    std::cout << "mirrored int64: "
              << *(volatile uint64_t*)virt_addr2 << std::endl;
    *(volatile uint64_t*)virt_addr = 1;
    msync(virt_addr, page_size, MS_SYNC);
    msync(virt_addr2, page_size, MS_INVALIDATE);
    std::cout << "mirrored int64: "
              << *(volatile uint64_t*)virt_addr2 << std::endl;

    test_write(virt_addr, 1 << 25);
    test_write(virt_addr2, 1 << 25);

    return 0;
}
