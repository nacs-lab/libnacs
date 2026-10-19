// Equivalence harness: device output of one serialized sequence, through ANY build of
// libnacs-utils.dll (loaded by explicit path, plain C API only, so the same program runs
// against the live DLL backup and a new build).
//
//   nacs-seq-devdump <libnacs-utils.dll> <seq.bin> <config.yml> <globals.txt|-> <outdir>
//
// Writes into outdir: builder.txt, seq.txt, seq_opt.txt (compiled IR dumps), and per
// basic sequence k (in run order): bseq<k>_zynq.bin (FPGA1 bytecode), bseq<k>_nidaq.bin
// (NiDAQ samples), bseq<k>_len.txt; msgs_*.bin = the manager's messages after each call.
//
// No hardware: the config must give the AWG devices backend `dummy` and the FPGA a dead
// local url; start() is never called (the FPGA sends its bytecode only in start()).
// Build: gcc -O2 -o nacs-seq-devdump.exe nacs_seq_devdump.c   (msys2 MINGW64)

#include <windows.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef void *(*new_manager_t)(void);
typedef void (*load_config_t)(void*, const char*);
typedef void (*enable_t)(void*, int);
typedef void *(*create_seq_t)(void*, const uint8_t*, size_t);
typedef uint8_t *(*take_msgs_t)(void*, size_t*);
typedef void (*seq_void_t)(void*);
typedef uint32_t (*post_run_t)(void*);
typedef void (*set_global_t)(void*, uint32_t, double);
typedef uint64_t (*len_t)(void*);
typedef const uint8_t *(*dump_t)(void*, size_t*);
typedef const uint8_t *(*zynq_bc_t)(void*, const char*, size_t*);
typedef const double *(*nidaq_t)(void*, const char*, size_t*);

static HMODULE dll;
static void *sym(const char *name)
{
    void *p = (void*)GetProcAddress(dll, name);
    if (!p) {
        fprintf(stderr, "missing symbol %s\n", name);
        exit(3);
    }
    return p;
}

static void write_file(const char *dir, const char *name, const void *data, size_t sz)
{
    char path[2048];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "cannot write %s\n", path);
        exit(4);
    }
    if (sz)
        fwrite(data, 1, sz, f);
    fclose(f);
}

static void *mgr;
static take_msgs_t take_msgs;
static int nmsg = 0;
static void save_msgs(const char *dir, const char *what)
{
    size_t sz = 0;
    uint8_t *m = take_msgs(mgr, &sz);
    char name[256];
    snprintf(name, sizeof(name), "msgs_%02d_%s.bin", nmsg++, what);
    write_file(dir, name, m, sz);
}

int main(int argc, char **argv)
{
    if (argc < 6) {
        fprintf(stderr, "usage: nacs-seq-devdump <dll> <seq.bin> <config.yml> <globals|-> <outdir>\n");
        return 1;
    }
    const char *out = argv[5];
    CreateDirectoryA(out, NULL);
    dll = LoadLibraryExA(argv[1], NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
    if (!dll) {
        fprintf(stderr, "cannot load %s (error %lu)\n", argv[1], GetLastError());
        return 2;
    }
    char loaded[2048];
    GetModuleFileNameA(dll, loaded, sizeof(loaded));
    printf("loaded %s\n", loaded);

    new_manager_t new_manager = (new_manager_t)sym("nacs_seq_new_manager");
    load_config_t load_config = (load_config_t)sym("nacs_seq_manager_load_config_file");
    enable_t enable_dump = (enable_t)sym("nacs_seq_manager_enable_dump");
    create_seq_t create_seq = (create_seq_t)sym("nacs_seq_manager_create_sequence");
    take_msgs = (take_msgs_t)sym("nacs_seq_manager_take_messages");
    seq_void_t init_run = (seq_void_t)sym("nacs_seq_manager_expseq_init_run");
    seq_void_t pre_run = (seq_void_t)sym("nacs_seq_manager_expseq_pre_run");
    post_run_t post_run = (post_run_t)sym("nacs_seq_manager_expseq_post_run");
    set_global_t set_global = (set_global_t)sym("nacs_seq_manager_expseq_set_global");
    len_t cur_len = (len_t)sym("nacs_seq_manager_expseq_cur_bseq_length");
    dump_t builder_dump = (dump_t)sym("nacs_seq_manager_expseq_get_builder_dump");
    dump_t seq_dump = (dump_t)sym("nacs_seq_manager_expseq_get_seq_dump");
    dump_t opt_dump = (dump_t)sym("nacs_seq_manager_expseq_get_seq_opt_dump");
    zynq_bc_t zynq_bc = (zynq_bc_t)sym("nacs_seq_manager_expseq_get_zynq_bytecode");
    nidaq_t nidaq = (nidaq_t)sym("nacs_seq_manager_expseq_get_nidaq_data");

    mgr = new_manager();
    load_config(mgr, argv[3]);
    save_msgs(out, "config");
    enable_dump(mgr, 1);

    FILE *f = fopen(argv[2], "rb");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", argv[2]);
        return 1;
    }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t *buf = (uint8_t*)malloc(n);
    if (fread(buf, 1, n, f) != (size_t)n)
        return 1;
    fclose(f);

    void *seq = create_seq(mgr, buf, (size_t)n);
    save_msgs(out, "create");
    if (!seq) {
        fprintf(stderr, "create_sequence failed (see msgs)\n");
        return 5;
    }
    size_t sz;
    const uint8_t *d;
    d = builder_dump(seq, &sz);
    write_file(out, "builder.txt", d, d ? sz : 0);
    d = seq_dump(seq, &sz);
    write_file(out, "seq.txt", d, d ? sz : 0);
    d = opt_dump(seq, &sz);
    write_file(out, "seq_opt.txt", d, d ? sz : 0);

    if (strcmp(argv[4], "-") != 0) {
        FILE *g = fopen(argv[4], "r");
        unsigned id;
        double val;
        while (g && fscanf(g, "%u %lf", &id, &val) == 2)
            set_global(seq, id, val);
        if (g)
            fclose(g);
    }
    init_run(seq);
    save_msgs(out, "init");
    int k = 0;
    while (k < 64) {
        pre_run(seq);
        char what[64];
        snprintf(what, sizeof(what), "pre_run%d", k);
        save_msgs(out, what);
        char name[256];
        d = zynq_bc(seq, "FPGA1", &sz);
        snprintf(name, sizeof(name), "bseq%d_zynq.bin", k);
        write_file(out, name, d, d ? sz : 0);
        const double *nd = nidaq(seq, "NiDAQ", &sz);
        snprintf(name, sizeof(name), "bseq%d_nidaq.bin", k);
        write_file(out, name, nd, nd ? sz * sizeof(double) : 0);
        char len[64];
        int ln = snprintf(len, sizeof(len), "%llu\n", (unsigned long long)cur_len(seq));
        snprintf(name, sizeof(name), "bseq%d_len.txt", k);
        write_file(out, name, len, ln);
        uint32_t next = post_run(seq);
        snprintf(what, sizeof(what), "post_run%d", k);
        save_msgs(out, what);
        k++;
        if (next == 0)
            break;
    }
    printf("ok: %d basic sequences\n", k);
    fflush(stdout);
    // Skip DLL teardown (ZeroMQ's exit-time cleanup asserts once Winsock is gone).
    TerminateProcess(GetCurrentProcess(), 0);
    return 0;
}
