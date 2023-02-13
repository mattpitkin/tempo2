#include <stdio.h>
#include <assert.h>

int main(int argc, char** argv) {
    assert(argc == 3);
    char* fn = argv[1];
    char* arrn = argv[2];
    FILE* f = fopen(fn, "rb");
    printf("const char %s[] = {\n",arrn);
    unsigned long n = 0;
    while(!feof(f)) {
        unsigned char c;
        if(fread(&c, 1, 1, f) == 0) break;
        printf("0x%.2X,", (int)c);
        ++n;
        if(n % 10 == 0) printf("\n");
    }
    fclose(f);
    printf("0x00};\n");
    return 0;
}
