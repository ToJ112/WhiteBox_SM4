#include "WBMatrix.h"
#include "sm4.h"

#define GET32(pc)  (\
((uint32_t)(pc)[0] << 24) ^\
((uint32_t)(pc)[1] << 16) ^\
((uint32_t)(pc)[2] <<  8) ^\
((uint32_t)(pc)[3]))

#define PUT32(st, ct)\
(ct)[0] = (uint8_t)((st) >> 24);\
(ct)[1] = (uint8_t)((st) >> 16);\
(ct)[2] = (uint8_t)((st) >>  8);\
(ct)[3] = (uint8_t)(st)

M32 L_matrix = {
        .M[0] = 0xA0202080,
        .M[1] = 0x50101040,
        .M[2] = 0x28080820,
        .M[3] = 0x14040410,
        .M[4] = 0xA020208,
        .M[5] = 0x5010104,
        .M[6] = 0x2808082,
        .M[7] = 0x1404041,
        .M[8] = 0x80A02020,
        .M[9] = 0x40501010,
        .M[10] = 0x20280808,
        .M[11] = 0x10140404,
        .M[12] = 0x80A0202,
        .M[13] = 0x4050101,
        .M[14] = 0x82028080,
        .M[15] = 0x41014040,
        .M[16] = 0x2080A020,
        .M[17] = 0x10405010,
        .M[18] = 0x8202808,
        .M[19] = 0x4101404,
        .M[20] = 0x2080A02,
        .M[21] = 0x1040501,
        .M[22] = 0x80820280,
        .M[23] = 0x40410140,
        .M[24] = 0x202080A0,
        .M[25] = 0x10104050,
        .M[26] = 0x8082028,
        .M[27] = 0x4041014,
        .M[28] = 0x202080A,
        .M[29] = 0x1010405,
        .M[30] = 0x80808202,
        .M[31] = 0x40404101
};



Aff32 M[32][3];
Aff32 C[32];
Aff32 D[32];
Aff32 SE[4];
Aff32 FE[4];
uint32_t Table[32][4][256];
uint32_t Table_part2[32][4][256];  //32轮 part2有四个这种表，8bit输入32bit输出   复合了嵌入rk，sbox，L移位
uint32_t Table_SL[32][4][256];     //经过sbox和L循环
uint8_t Table_addIn_part2[32][4][256];
uint32_t Mask[32][4];              //32轮，每轮4个32bit mask
Aff8 Eij[32][4];
Aff8 Eij_inv[32][4];
Aff32 Ei[32];
Aff32 MB[32][4];
Aff32 MB_inv[32][4];
uint8_t Out_part2[32][4][8][16];          //32轮，每轮4个表，暂时定4*(32/4)=32个out一轮要用，out编码输入是4bit即16种可能值
uint8_t In_part2[32][4][8][16];          //32轮，每轮4个表，暂时定4*(32/4)=32个in一轮要用，in编码输入是4bit即16种可能值
uint8_t Out_part1[32][8][16];
uint8_t In_part1[32][8][16];

void printstate(unsigned char* in);
void wbsm4_gen(uint8_t* key);
void wbsm4_encrypt(unsigned char IN[], unsigned char OUT[]);
