#include <nacs-utils/fd_utils.h>
#include <nacs-utils/mem.h>
#include <nacs-utils/timer.h>

#include <unistd.h>
#include <sys/mman.h>
#include <iostream>
#include <iomanip>

using namespace NaCs;

void
test_write(volatile void *ptr, uint32_t nrun)
{
    tic();
    for (uint32_t i = 0;i < nrun;i++) {
        mem_write64(ptr, i);
    }
    auto time = toc();
    std::cout << "Time per write: " << std::setprecision(4)
              << double(time) / double(nrun) / 1e3
              << " us" << std::endl;
}

int
main()
{
    uint32_t page_size = getpagesize();

    void *virt_addr = mmap(nullptr, page_size, PROT_READ | PROT_WRITE,
                           MAP_SHARED | MAP_ANONYMOUS | MAP_LOCKED |
                           MAP_POPULATE, -1, 0);
    void *phy_addr = getPhyAddr(virt_addr);

    void *virt_addr2 = mapFile("/dev/mem", off_t(phy_addr), page_size);
    strcpy((char*)virt_addr2, "random string");
    printf("%s\n", virt_addr);
    strcpy((char*)virt_addr2, "random string22222");
    printf("%s\n", virt_addr);

    // strcpy((char*)virt_addr, "random string");
    // printf("%s\n", virt_addr2);
    // strcpy((char*)virt_addr, "random string22222");
    // printf("%s\n", virt_addr2);

    std::cout << "page_size: 0x" << std::hex << page_size << std::endl;
    std::cout << "phy_addr : 0x" << std::hex << phy_addr << std::endl;
    std::cout << "virt_addr: 0x" << std::hex << virt_addr << std::endl;
    std::cout << "virt_addr2: 0x" << std::hex << virt_addr2 << std::endl;
    std::cout << "phy_addr2: 0x" << std::hex
              << getPhyAddr(virt_addr2) << std::endl;

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
