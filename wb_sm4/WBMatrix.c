#include "WBMatrix.h"

unsigned int randseed;
uint8_t idM8[8] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
uint32_t idM32[32] = { 0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000, 0x2000000, 0x1000000, 0x800000, 0x400000, 0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1 };
//8bit internal xor table
int xor [] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0,
1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0,
1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0,
1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0,
1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0,
1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1,
0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0 };

#ifndef GET_ULONG_BE
#define GET_ULONG_BE(n,b,i)                             \
{                                                       \
    (n) = ( (unsigned long) (b)[(i)    ] << 24 )        \
        | ( (unsigned long) (b)[(i) + 1] << 16 )        \
        | ( (unsigned long) (b)[(i) + 2] <<  8 )        \
        | ( (unsigned long) (b)[(i) + 3]       );       \
}
#endif
void SetRandSeed(unsigned int seed)
{
    randseed = seed;
}
int xorU8(uint8_t n)// uint8_t internal xor
{
    if (xor [n]) return 1;
    else return 0;
}
int xorU16(uint16_t n)// uint16_t internal xor
{
    uint8_t temp = 0;
    uint8_t* u = (uint8_t*)&n;
    temp = (*u) ^ (*(u + 1));
    if (xorU8(temp)) return 1;
    else return 0;
}
int xorU32(uint32_t n)// uint32_t internal xor
{
    uint16_t temp = 0;
    uint16_t* u = (uint16_t*)&n;
    temp = (*u) ^ (*(u + 1));
    if (xorU16(temp)) return 1;
    else return 0;
}
/*
void V32toV4(uint32_t n,uint8_t* b[8]) {
    for (int i = 0; i < 8; i++) {
        b[i] = (unsigned char)((n >> (28 - i * 4))&0x0F);
    }
}
void V4toV32(uint32_t* n, uint8_t b[8]) {
    n = ((uint32_t) (b)[0]<<28)   \
      | ((uint32_t)(b)[1] << 24)  \
      | ((uint32_t)(b)[2] << 20)  \
      | ((uint32_t)(b)[3] << 16)  \
      | ((uint32_t)(b)[4] << 12)  \
      | ((uint32_t)(b)[5] << 8)   \
      | ((uint32_t)(b)[6] << 44)  \
      | ((uint32_t)(b)[7]      )  \
}*/
void initV8(V8* Vec)//initial Vector 8*1
{
    (*Vec).V = 0;
}
void initV32(V32* Vec)//initial Vector 32*1
{
    (*Vec).V = 0;
}
void MatMulVecM32(M32 Mat, V32 Vec, V32* ans)//matrix * vector -> vector 32*1
{
    int i;
    initV32(ans);
    for (i = 0; i < 32; i++)
    {
        if (xorU32(Mat.M[i] & Vec.V)) (*ans).V ^= idM32[i];
    }
}
void MatMulVecM8(M8 Mat, V8 Vec, V8* ans)//matrix * vector -> vector 8*1
{
    int i;
    initV8(ans);
    for (i = 0; i < 8; i++)
    {
        if (xorU8(Mat.M[i] & Vec.V)) (*ans).V ^= idM8[i];
    }
}
uint32_t MatMulNumM32(M32 Mat, uint32_t n)//matrix * number -> number 32bits
{
    int i;
    uint32_t temp = 0;
    for (i = 0; i < 32; i++)
    {
        if (xorU32(Mat.M[i] & n)) temp ^= idM32[i];
    }
    return temp;
}
uint32_t VecAddVecV32(uint32_t Vec1, uint32_t Vec2) {    //32bitÏòÁ¿Ïà³Ë
    return (Vec1 ^ Vec2);
}
void identityM8(M8* Mat)//identity matrix 8*8
{
    int i;
    for (i = 0; i < 8; i++)
    {
        (*Mat).M[i] = idM8[i];
    }
}
void identityM32(M32* Mat)//identity matrix 32*32
{
    int i;
    for (i = 0; i < 32; i++)
    {
        (*Mat).M[i] = idM32[i];
    }
}
void randV8(V8* Vec)//randomize Vector 8*1
{
    InitRandom((randseed++) ^ (unsigned int)time(NULL));
    (*Vec).V = cus_random();
}
void randV32(V32* Vec)//randomize Vector 32*1
{
    uint16_t* v = (uint16_t*)&((*Vec).V);
    InitRandom((randseed++) ^ (unsigned int)time(NULL));
    *(v + 1) = cus_random();
    *v = cus_random();
}
void randM8(M8* Mat)//randomize Matrix 8*8 
{
    int i;
    InitRandom((randseed++) ^ ((unsigned int)time(NULL)));
    for (i = 0; i < 8; i++)
    {
        (*Mat).M[i] = cus_random();
    }
}
void randM32(M32* Mat)//randomize Matrix 32*32 
{
    int i;
    InitRandom((randseed++) ^ ((unsigned int)time(NULL)));
    for (i = 0; i < 32; i++)
    {
        (*Mat).M[i] = cus_random();
    }
}
void copyM8(M8 Mat1, M8* Mat2)
{
    int i;
    for (i = 0; i < 8; i++)
    {
        (*Mat2).M[i] = Mat1.M[i];
    }
}
void copyM32(M32 Mat1, M32* Mat2)
{
    int i;
    for (i = 0; i < 32; i++)
    {
        (*Mat2).M[i] = Mat1.M[i];
    }
}
int isequalM32(M32 Mat1, M32 Mat2)
{
    int i;
    int flag = 1;
    for (i = 0; i < 32; i++)
    {
        if (Mat1.M[i] != Mat2.M[i])
        {
            flag = 0;
            break;
        }
    }
    return flag;
}
int isinvertM32(M32 Mat)//Invertible Matrix?
{
    int i, j, k;
    uint32_t temp;
    int flag;
    for (i = 0; i < 32; i++)
    {
        if ((Mat.M[i] & idM32[i]) == idM32[i])
        {
            for (j = i + 1; j < 32; j++)
            {
                if ((Mat.M[j] & idM32[i]) == idM32[i])
                {
                    Mat.M[j] ^= Mat.M[i];
                }
            }
        }
        else
        {
            flag = 1;
            for (j = i + 1; j < 32; j++)
            {
                if ((Mat.M[j] & idM32[i]) == idM32[i])
                {
                    temp = Mat.M[i];
                    Mat.M[i] = Mat.M[j];
                    Mat.M[j] = temp;
                    flag = 0;
                    break;
                }
            }
            if (flag) return 0;
            for (k = i + 1; k < 32; k++)
            {
                if ((Mat.M[k] & idM32[i]) == idM32[i])
                {
                    Mat.M[k] ^= Mat.M[i];
                }
            }
        }
    }
    if (Mat.M[31] == idM32[31]) return 1;
    else return 0;
}

uint32_t affineU32(Aff32 aff, uint32_t arr)//32bits affine transformation
{
    V32 mul_vec, ans_vec;
    initV32(&ans_vec);
    mul_vec.V = arr;
    MatMulVecM32(aff.Mat, mul_vec, &ans_vec);//mul
    return ans_vec.V ^ aff.Vec.V;//add
}
void genMatpairM8(M8* Mat, M8* Mat_inv)//generate 8*8 invertible matrix and its inverse matrix
{
    int i, j, t, k;
    int p, q;
    M8 tempMat;
    M8 resultMat;
    uint8_t temp;
    uint8_t trail[64][3];// generate trail
    int flag = 0;
    int times = 0;
    int invertible = 1;
    InitRandom((randseed++) ^ ((unsigned int)time(NULL)));
    identityM8(Mat);
    identityM8(Mat_inv);
    randM8(&tempMat);
    copyM8(tempMat, &resultMat);
    for (i = 0; i < 8; i++)//diagonal = 1?
    {
        if ((tempMat.M[i] & idM8[i]) == idM8[i])
        {
            for (j = i + 1; j < 8; j++)
            {
                if ((tempMat.M[j] & idM8[i]) == idM8[i])
                {
                    tempMat.M[j] ^= tempMat.M[i];

                    (*Mat_inv).M[j] ^= (*Mat_inv).M[i];

                    trail[times][0] = 1;
                    trail[times][1] = j;
                    trail[times][2] = i;
                    times++;
                }
            }
        }
        else// swap to find 1
        {
            flag = 1;
            for (j = i + 1; j < 8; j++)
            {
                if ((tempMat.M[j] & idM8[i]) == idM8[i])
                {
                    temp = tempMat.M[i];
                    tempMat.M[i] = tempMat.M[j];
                    tempMat.M[j] = temp;

                    flag = 0;

                    temp = (*Mat_inv).M[i];
                    (*Mat_inv).M[i] = (*Mat_inv).M[j];
                    (*Mat_inv).M[j] = temp;

                    trail[times][0] = 0;
                    trail[times][1] = j;
                    trail[times][2] = i;
                    times++;
                    break;
                }
            }
            if (flag) //can not find 1 which means not invertible
            {
                invertible = 0;
                if (i < 7)
                {
                    p = i + 1 + cus_random() % (7 - i);//swap
                    temp = tempMat.M[p];
                    tempMat.M[p] = tempMat.M[i];
                    tempMat.M[i] = temp;
                    temp = (*Mat_inv).M[p];
                    (*Mat_inv).M[p] = (*Mat_inv).M[i];
                    (*Mat_inv).M[i] = temp;
                    trail[times][0] = 0;
                    trail[times][1] = p;
                    trail[times][2] = i;
                    times++;
                    for (t = i + 1; t < 8; t++)
                    {
                        if (cus_random() % 2)
                        {
                            tempMat.M[t] ^= tempMat.M[i];
                            (*Mat_inv).M[t] ^= (*Mat_inv).M[i];
                            trail[times][0] = 1;
                            trail[times][1] = t;
                            trail[times][2] = i;
                            times++;
                        }
                    }
                }
            }
            else //can still contiune
            {
                for (k = i + 1; k < 8; k++)
                {
                    if ((tempMat.M[k] & idM8[i]) == idM8[i])
                    {
                        tempMat.M[k] ^= tempMat.M[i];

                        (*Mat_inv).M[k] ^= (*Mat_inv).M[i];

                        trail[times][0] = 1;
                        trail[times][1] = k;
                        trail[times][2] = i;
                        times++;
                    }
                }
            }
        }
    }
    if (!invertible)//not invertible
    {
        for (t = 7; t >= 0; t--)
        {
            for (j = t - 1; j >= 0; j--)
            {
                if ((tempMat.M[j] & idM8[t]) == idM8[t])
                {
                    tempMat.M[j] ^= tempMat.M[t];
                    (*Mat_inv).M[j] ^= (*Mat_inv).M[t];
                    trail[times][0] = 1;
                    trail[times][1] = j;
                    trail[times][2] = t;
                    times++;
                }
            }
        }

        for (j = times - 1; j >= 0; j--)//generate inverse matrix
        {
            if (trail[j][0])//add
            {
                (*Mat).M[trail[j][1]] ^= (*Mat).M[trail[j][2]];
            }
            else//swap
            {
                temp = (*Mat).M[trail[j][1]];
                (*Mat).M[trail[j][1]] = (*Mat).M[trail[j][2]];
                (*Mat).M[trail[j][2]] = temp;
            }
        }
    }
    else//invertible 
    {
        for (i = 7; i >= 0; i--)
        {
            for (j = i - 1; j >= 0; j--)
            {
                if ((tempMat.M[j] & idM8[i]) == idM8[i])
                {
                    tempMat.M[j] ^= tempMat.M[i];

                    (*Mat_inv).M[j] ^= (*Mat_inv).M[i];
                }
            }
        }
        copyM8(resultMat, Mat);
    }
}
void genMatpairM32(M32* Mat, M32* Mat_inv)//generate 32*32 invertible matrix and its inverse matrix
{
    int i, j, t, k;
    int p, q;
    M32 tempMat;
    M32 resultMat;
    uint32_t temp;
    uint8_t trail[1024][3];// generate trail
    int flag = 0;
    int times = 0;
    int invertible = 1;
    InitRandom((randseed++) ^ ((unsigned int)time(NULL)));
    identityM32(Mat);
    identityM32(Mat_inv);
    randM32(&tempMat);
    copyM32(tempMat, &resultMat);
    for (i = 0; i < 32; i++)//diagonal = 1?
    {
        if ((tempMat.M[i] & idM32[i]) == idM32[i])
        {
            for (j = i + 1; j < 32; j++)
            {
                if ((tempMat.M[j] & idM32[i]) == idM32[i])
                {
                    tempMat.M[j] ^= tempMat.M[i];

                    (*Mat_inv).M[j] ^= (*Mat_inv).M[i];

                    trail[times][0] = 1;
                    trail[times][1] = j;
                    trail[times][2] = i;
                    times++;
                }
            }
        }
        else// swap to find 1
        {
            flag = 1;
            for (j = i + 1; j < 32; j++)
            {
                if ((tempMat.M[j] & idM32[i]) == idM32[i])
                {
                    temp = tempMat.M[i];
                    tempMat.M[i] = tempMat.M[j];
                    tempMat.M[j] = temp;

                    flag = 0;

                    temp = (*Mat_inv).M[i];
                    (*Mat_inv).M[i] = (*Mat_inv).M[j];
                    (*Mat_inv).M[j] = temp;

                    trail[times][0] = 0;
                    trail[times][1] = j;
                    trail[times][2] = i;
                    times++;
                    break;
                }
            }
            if (flag) //can not find 1 which means not invertible
            {
                invertible = 0;
                if (i < 31)
                {
                    p = i + 1 + cus_random() % (31 - i);//swap
                    temp = tempMat.M[p];
                    tempMat.M[p] = tempMat.M[i];
                    tempMat.M[i] = temp;
                    temp = (*Mat_inv).M[p];
                    (*Mat_inv).M[p] = (*Mat_inv).M[i];
                    (*Mat_inv).M[i] = temp;
                    trail[times][0] = 0;
                    trail[times][1] = p;
                    trail[times][2] = i;
                    times++;
                    for (t = i + 1; t < 32; t++)
                    {
                        if (cus_random() % 2)
                        {
                            tempMat.M[t] ^= tempMat.M[i];
                            (*Mat_inv).M[t] ^= (*Mat_inv).M[i];
                            trail[times][0] = 1;
                            trail[times][1] = t;
                            trail[times][2] = i;
                            times++;
                        }
                    }
                }
            }
            else //can still contiune
            {
                for (k = i + 1; k < 32; k++)
                {
                    if ((tempMat.M[k] & idM32[i]) == idM32[i])
                    {
                        tempMat.M[k] ^= tempMat.M[i];

                        (*Mat_inv).M[k] ^= (*Mat_inv).M[i];

                        trail[times][0] = 1;
                        trail[times][1] = k;
                        trail[times][2] = i;
                        times++;
                    }
                }
            }
        }
    }
    if (!invertible)//not invertible
    {
        for (t = 31; t >= 0; t--)
        {
            for (j = t - 1; j >= 0; j--)
            {
                if ((tempMat.M[j] & idM32[t]) == idM32[t])
                {
                    tempMat.M[j] ^= tempMat.M[t];
                    (*Mat_inv).M[j] ^= (*Mat_inv).M[t];
                    trail[times][0] = 1;
                    trail[times][1] = j;
                    trail[times][2] = t;
                    times++;
                }
            }
        }

        for (j = times - 1; j >= 0; j--)//generate inverse matrix
        {
            if (trail[j][0])//add
            {
                (*Mat).M[trail[j][1]] ^= (*Mat).M[trail[j][2]];
            }
            else//swap
            {
                temp = (*Mat).M[trail[j][1]];
                (*Mat).M[trail[j][1]] = (*Mat).M[trail[j][2]];
                (*Mat).M[trail[j][2]] = temp;
            }
        }
    }
    else//invertible 
    {
        for (i = 31; i >= 0; i--)
        {
            for (j = i - 1; j >= 0; j--)
            {
                if ((tempMat.M[j] & idM32[i]) == idM32[i])
                {
                    tempMat.M[j] ^= tempMat.M[i];

                    (*Mat_inv).M[j] ^= (*Mat_inv).M[i];
                }
            }
        }
        copyM32(resultMat, Mat);
    }
}
void genaffinepairM8(Aff8* aff, Aff8* aff_inv)//generate a pair of affine
{
    genMatpairM8(&(aff->Mat), &(aff_inv->Mat));
    randV8(&(aff->Vec));
    MatMulVecM8((*aff_inv).Mat, (*aff).Vec, &(aff_inv->Vec));
}
void genaffinepairM32(Aff32* aff, Aff32* aff_inv)//generate a pair of affine
{
    genMatpairM32(&(aff->Mat), &(aff_inv->Mat));
    randV32(&(aff->Vec));
    MatMulVecM32((*aff_inv).Mat, (*aff).Vec, &(aff_inv->Vec));
}
void MatrixcomM8to32(M8 m1, M8 m2, M8 m3, M8 m4, M32* mat)//diagonal matrix concatenation, four 8*8 -> 32*32
{
    int i;
    int j = 0;
    uint8_t* m;
    initM32(mat);
    for (i = 0; i < 8; i++)
    {
        m = (uint8_t*)&(*mat).M[j];
        *(m + 3) = m1.M[i];
        j++;
    }
    for (i = 0; i < 8; i++)
    {
        m = (uint8_t*)&(*mat).M[j];
        *(m + 2) = m2.M[i];
        j++;
    }
    for (i = 0; i < 8; i++)
    {
        m = (uint8_t*)&(*mat).M[j];
        *(m + 1) = m3.M[i];
        j++;
    }
    for (i = 0; i < 8; i++)
    {
        m = (uint8_t*)&(*mat).M[j];
        *m = m4.M[i];
        j++;
    }
}
void VectorcomV8to32(V8 v1, V8 v2, V8 v3, V8 v4, V32* vec)//4 vectors concatenation
{
    uint8_t* v;
    v = (uint8_t*)&(*vec).V;
    *(v + 3) = v1.V;
    *(v + 2) = v2.V;
    *(v + 1) = v3.V;
    *v = v4.V;
}
void affinecomM8to32(Aff8 aff1, Aff8 aff2, Aff8 aff3, Aff8 aff4, Aff32* aff)//diagonal affine concatenation, four 8*8 -> 32*32
{
    MatrixcomM8to32(aff1.Mat, aff2.Mat, aff3.Mat, aff4.Mat, &(aff->Mat));
    VectorcomV8to32(aff1.Vec, aff2.Vec, aff3.Vec, aff4.Vec, &(aff->Vec));
}