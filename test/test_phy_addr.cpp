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
                           MAP_SHARED | MAP_ANONYMOUS | MAP_LOCKED, -1, 0);
    auto virt_pfn = uintptr_t(virt_addr) / page_size;
    int fd = open("/proc/self/pagemap", O_RDONLY);
    uint64_t page_info;
    pread(fd, &page_info, sizeof(page_info), virt_pfn * sizeof(uint64_t));
    uint64_t page_info2;
    pread(fd, &page_info2, sizeof(page_info2), virt_pfn * sizeof(uint64_t));
    close(fd);
    uint64_t phy_pfn = page_info & ((1ll << 55) - 1);
    uint64_t phy_addr = phy_pfn * page_size;

    volatile void *virt_addr2 = mapFile("/dev/mem", off_t(phy_addr), page_size);
    strcpy((char*)virt_addr, "random string");
    printf("%s\n", virt_addr2);

    std::cout << "page_size: 0x" << std::hex << page_size << std::endl;
    std::cout << "phy_addr : 0x" << std::hex << phy_addr << std::endl;
    std::cout << "virt_addr: 0x" << std::hex << virt_addr << std::endl;
    std::cout << "page_info: 0x" << std::hex << page_info << std::endl;
    std::cout << "virt_addr2: 0x" << std::hex << (void*)virt_addr2 << std::endl;

    test_write(virt_addr, 1 << 25);
    test_write(virt_addr2, 1 << 25);

    return 0;
}
