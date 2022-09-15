#ifndef _HWBMATRIX_H_
#define _HWBMATRIX_H_

#include "structure.h"
#include "random.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>



void SetRandSeed(unsigned int seed);//Set random seed
int xorU32(uint32_t n);
void V32toV4(uint32_t n, uint8_t* b[8]);
void V4toV32(uint32_t* n, uint8_t b[8]);
void initV32(V32* Vec);
void MatMulVecM32(M32 Mat, V32 Vec, V32* ans);

uint32_t MatMulNumM32(M32 Mat, uint32_t n);   //32*32矩阵与32维向量相乘
uint32_t VecAddVecV32(uint32_t Vec1, uint32_t Vec2);   //32bit向量的异或
void identityM32(M32* Mat);
void randV32(V32* Vec);
void randM32(M32* Mat);
void copyM32(M32 Mat1, M32* Mat2);
int isequalM32(M32 Mat1, M32 Mat2);
int isinvertM32(M32 Mat);
uint32_t affineU32(Aff32 aff, uint32_t arr);
void genMatpairM32(M32* Mat, M32* Mat_inv);
void genaffinepairM32(Aff32* aff, Aff32* aff_inv);


#endif