#include "wbsm4.h"
#include "sbox.h"

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


uint32_t Table_part2[32][4][256];  //32�� part2���ĸ����ֱ�8bit����32bit���   ������Ƕ��rk��sbox��L��λ
uint32_t Table_SL[32][4][256];     //����sbox��Lѭ��
uint8_t Table_addIn_part2[32][4][256];
uint32_t Mask[32][4];              //32�֣�ÿ��4��32bit mask
Aff8 Eij[32][4];
Aff8 Eij_inv[32][4];
Aff32 Ei[32];
Aff32 MB[32][4];
Aff32 MB_inv[32][4];
uint8_t Out_part2[32][4][8][16];          //32�֣�ÿ��4������ʱ��4*(32/4)=32��outһ��Ҫ�ã�out����������4bit��16�ֿ���ֵ
uint8_t In_part2[32][4][8][16];          //32�֣�ÿ��4������ʱ��4*(32/4)=32��inһ��Ҫ�ã�in����������4bit��16�ֿ���ֵ
uint8_t Out_part1[32][8][16];
uint8_t In_part1[32][8][16];

void printstate(unsigned char* in)
{
    int i;
    for (i = 0; i < 16; i++)
    {
        printf("%.2X", in[i]);
    }
    printf("\n");
}
void randomOutIn(unsigned char Out[32][4][8][16], unsigned char In[32][4][8][16]) {//ֻ��һ����Ϊ��

    //��ʼ��
    srand((unsigned int)time(NULL));		//ʱ�䲥��
    for (int r = 0; r < 32; r++)					//ǰ9��
        for (int i = 0; i < 4; i++) {			//ÿһ����4���л��
            for (int j = 0; j < 8; j++) {		//ÿһ���л�ϣ����õ���Out����
                for (int k = 0; k < 16; k++) {	//��СΪ16��һά��������4bit��4bit���������
                    Out[r][i][j][k] = k;

                }
            }
        }

    //����û�����Out 
    for (int r = 0; r < 32; r++)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                for (int k = 0; k < 16; k++) {
                    int randNum = rand() % 16;
                	
                    int t = Out[r][i][j][k];
                    Out[r][i][j][k] = Out[r][i][j][randNum];
                    Out[r][i][j][randNum] = t;
                }
            }
        }

    //����In 
    for (int r = 0; r < 32; r++)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                for (int k = 0; k < 16; k++) {
                    In[r][i][j][Out[r][i][j][k]] = k;
                }
            }
        }
    return;
}


void wbsm4_gen_part1Table() {

}

void wbsm4_gen_part2Table(uint8_t* key)
 {
    sm4_context ctx;
    sm4_setkey_enc(&ctx, key);   //��Կ��չ
    InitRandom(((unsigned int)time(NULL)));
    //����32�֣�ÿ��4����  �����������
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            Mask[i][j] = cus_random();
        }
    }
    for (int i = 0; i < 32; i++)
    {
        //affine E 
        for (int j = 0; j < 4; j++)
        {
            genaffinepairM8(&Eij[i][j], &Eij_inv[i][j]);
        }

        // combine 4 E8 to 1 E32
        affinecomM8to32(Eij[i][0], Eij[i][1], Eij[i][2], Eij[i][3], &Ei[i]);
    }
    //����32�֣�ÿ���ĸ���  �������Ա������ͳ����Լ��������
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            genaffinepairM32(&MB[i][j], &MB_inv[i][j]);
        }
    }
    //����out�Լ�in�ķ����Ա���
    randomOutIn(Out_part2, In_part2);
    for (int i = 0; i < 4; i++)                         //��һ�ֽ���part2�������룬����������x���Ҳ��x
        for (int x = 0; x < 256; x++)
            Table_addIn_part2[0][i][x] = x;

    for (int r = 1; r < 32; r++) {                                //����s��֮ǰ��ȥ��part1�ķ����Ա���
        for(int i = 0; i < 4; i++){
            for (int x = 0; x < 256; x++) {
                uint8_t y = x;
                Table_addIn_part2[r][i][x] =
                    (In_part1[r][2 * i + 0][y >> 4 & 0xf] << 4) |
                    (In_part1[r][2 * i + 1][y >> 0 & 0xf] << 0);

            }
        }
    }


    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
                uint8_t temp_u8 = SBOX[Table_addIn_part2[i][j][x] ^ (ctx.sk[i] >> (24 - j * 8)) && 0xFF];   //���rk��sbox
                uint32_t temp_u32 = temp_u8 << (24 - j * 8);                      //8bit��չ��32bit
                Table_SL[i][j][x] = MatMulNumM32(L_matrix, temp_u32);             //���rk��sbox��L��λ
            }
        }
    }

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
               
                uint32_t temp_u32 = VecAddVecV32(Table_SL[i][j][x], Mask[i][j]);      //���rk��sbox��L��λ,+mask
                temp_u32 = MatMulNumM32(Ei[i].Mat, temp_u32);    //���rk��sbox��L��λ,+mask,+���Ծ���MB.Mat
                Table_part2[i][j][x] = VecAddVecV32(temp_u32, Ei[i].Vec.V);    //���rk��sbox��L��λ,+mask,+���Ծ���MB.Mat,+MB.Vec

            }
        }
    }
    for (int r = 0; r < 32; r++) {                                ////���rk��sbox��L��λ,+mask,+���Ծ���MB.Mat,+MB.Vec,+�����Ա���Out
        for (int i = 0; i < 4; i++) {
            for (int x = 0; x < 256; x++) {
                uint8_t y = x;
                uint32_t a = Table_part2[r][i][y];
                Table_part2[r][i][x] = 
                    (Out_part2[r][i][0][a >> 28 & 0xf] << 28) |
                    (Out_part2[r][i][1][a >> 24 & 0xf] << 24) |
                    (Out_part2[r][i][2][a >> 20 & 0xf] << 20) |
                    (Out_part2[r][i][3][a >> 16 & 0xf] << 16) |
                    (Out_part2[r][i][4][a >> 12 & 0xf] << 12) |
                    (Out_part2[r][i][5][a >> 8 & 0xf] << 8)   |
                    (Out_part2[r][i][6][a >> 4 & 0xf] << 4)   |
                    (Out_part2[r][i][7][a >> 0 & 0xf] << 0);

            }
        }
    }

}