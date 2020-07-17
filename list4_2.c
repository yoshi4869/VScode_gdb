#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define PI acos(-1) //3.1459
#define L 1000//number of transmit bits
#define SNR_STEP 60
/*データ型の定義*/
typedef struct
{
    double I;
    double Q;
}COMPLEX;
typedef struct
{
    /* data */
    COMPLEX *sig;
    unsigned int len;
}SEQ_SIG;
typedef struct
{
    /* data */
    unsigned short int *dat;
    unsigned int len;
}SEQ_DATA;

/*一様乱数生成*/
double _randU(void){
    return((double)rand()/(double)RAND_MAX);
}
/* ガウス雑音生成 */
/* 平均0，分散1のガウス雑音を一つ返す */
double _randN(void){
    double s, r, t;
    s = _randU();
    if(s == 0.0) s = 0.000000001;
    r = sqrt(-2.0*log(s));
    t = 2.0*PI*_randU();
    return(r*sin(t));
}
/*ランダムデータ発生*/
/* SEQ_DATA型変数dのメンバd.lenを長さとする0,1の2値系列を生成し
   アドレスd.datを先頭とするメモリ領域に保存する */
void _randData(SEQ_DATA d){
    unsigned int i;
    double x;
    for(i=0; i<d.len;i++){
        x = _randU();
        if(x >= 0.5)
            *(d.dat+i) = 1.0;
        else
            *(d.dat+i) = 0.0;
    }
}
/*AWGN発生*/
/* 長さがSEQ_SIG型変数nのメンバn.lenであり，電力Pnの複素AWGMを生成し，
   アドレスn.sigを先頭とするメモリ領域に保存する*/
void _awgn(SEQ_SIG n, double Pn){
    unsigned int i;
    for(i=0; i<n.len;i++){
        (n.sig+i)->I = _randN() * sqrt(Pn/2);
        (n.sig+i)->Q = _randN() * sqrt(Pn/2);
    }
}
/*雑音電力計算*/
double _SNRdB2noisePower(double c_dB){
    return pow(10,(-1)*c_dB/10);
}

/*BPSK変調器*/
/*SEQ_DATA型変数dataのメンバであるポインタdata.datを先頭アドレスとする
  メモリ領域に保存された長さdata.lenの2値データ系列にBPSK変調を施し
   結果をSEQ_SIG型変数sのメンバであるポインタs.sigを先頭アドレスとするメモリ領域に保存する*/
double _bpskMod(SEQ_SIG s, SEQ_DATA data){
    unsigned int i;
    if(s.len != data.len){
        printf("Error: _bepskMod, 長さが一致しません. \n");
        exit(1);
    }
    else{
        for(i=0;i<s.len;i++)
            (s.sig+i)->I = (*(data.dat+i) - 0.5)*2.0;
    }
}
/*ベクトル加算*/
/*2つの複素数列があり，SEQ_SIG型変数bとcで表される．
  つまり，長さはそれぞれb.lenおよび，c.lenであり，それぞれ先頭アドレスを
  b.sigおよび，c.sigとするメモリ領域に保存されている．
  これら2つの数列の要素ごとの足し算を計算し，得られた結果を先頭アドレスa.sigのメモリ領域に保存する*/
void _vectorSum(SEQ_SIG a, SEQ_SIG b, SEQ_SIG c){
    unsigned int i;
    if(a.len != b.len || b.len != c.len){
        printf("Error: _vectorSum, 長さが一致しません.\n");
        exit(1);
    }
    else{
        for(i=0;i<a.len;i++){
            (a.sig+i)->I = (b.sig+i)->I + (c.sig+i)->I;
            (a.sig+i)->Q = (b.sig+i)->Q + (c.sig+i)->Q;
        }
    }
}
/*BPSK復調器*/
/*SEQ_SIG型変数rsのメンバrs.lenを長さとする複素数列が先頭アドレスを
  rx.sigとするメモリ領域に保存されている．これにBPSK復調を施して
  得られた結果を先頭アドレスrdata.sigとするメモリ領域に保存する*/
void _bpskDem(SEQ_DATA rData, SEQ_SIG rs){
    unsigned int i;
    for(i=0;i<rs.len;i++){
        if((rs.sig+i)->I>0)
            *(rData.dat+i) = 1;
        else
            *(rData.dat+i) = 0;
    }
}
/*BER計算*/
/*SEQ_DATA型変数dataとrdataで表される2つのデータ系列がある．
  つまり，2つの系列の長さはそれぞれ，data.lenおよびrdata.lenであり，
  それぞれ先頭アドレスをdata.sigおよびrdata.sigとするメモリ領域に保存されている
  これら，2つのデータ系列を比較してBERを計算し，その値を返す*/
double _BER(SEQ_DATA data, SEQ_DATA rData){
    unsigned int i, sum = 0;
    double BER;
    if(data.len != rData.len){
        printf("Error: _BER() 長さが一致しません.\n");
        exit(1);
    }
    else{
        for(i=0;i<data.len;i++)
            sum += abs(*(data.dat+i) - *(rData.dat+i));
        BER = (double)sum/data.len;
    }
    return(BER);
}
/*得られたBERのファイルへ書き出し*/
/**SNRを先頭アドレスとするメモリ領域に保存されたsnr_step個のSN比と，
 * それに対応して*BERを先頭アドレスとするメモリ領域に保存されたBERの値を画面上及びファイルに書き出す*/
void _BERprint(unsigned int snr_step, double *SNRdB, double *BER){
    unsigned int i;
    FILE *fp;
    if((fp = fopen("BER.txt","w"))==NULL){
        printf("Error:Can not open your file \n");
        exit(1);
    }
    else{
        for(i=0;i<snr_step;i++){
            printf("%f [dB] BER ~ %f\n", *(SNRdB+i), *(BER+i));
            fprintf(fp,"%f %f\n",*(SNRdB+i),*(BER+i));
        }
        fclose(fp);
    }
}

int main(void){
    /*gnuplot*/
    FILE *gp,*gp2;
    /*複素数系列の型SEQ_SIGの定義*/
    SEQ_SIG txSymbol,rxSymbol, noise;
    SEQ_DATA txData,rxData;
    unsigned int i, x[L];
    double BER[SNR_STEP], SNRdB[SNR_STEP];

    txData.len = rxData.len = txSymbol.len = rxSymbol.len = noise.len = L;
    txData.dat = (unsigned short int *)malloc(sizeof(unsigned short int)*L);
    rxData.dat = (unsigned short int *)malloc(sizeof(unsigned short int)*L);
    txSymbol.sig = (COMPLEX *)malloc(sizeof(COMPLEX)*L);
    rxSymbol.sig = (COMPLEX *)malloc(sizeof(COMPLEX)*L);
    noise.sig = (COMPLEX *)malloc(sizeof(COMPLEX)*L);

    if(txData.dat == NULL||rxData.dat==NULL||txSymbol.sig==NULL||rxSymbol.sig==NULL||noise.sig==NULL){
        printf("Error:メモリを確保できません．\n");
        exit(1);
    }
    else{
            for(i=0; i<L; i++) x[i] = i;
            for(i=0; i<SNR_STEP; i++) SNRdB[i] = (double)i;
            for(i=0; i<SNR_STEP; i++){
            srand(time(NULL));

            _randData(txData); 

            _bpskMod(txSymbol,txData);

            _awgn(noise,_SNRdB2noisePower(SNRdB[i]));

            _vectorSum(rxSymbol,txSymbol,noise);

            _bpskDem(rxData,rxSymbol);

            BER[i]=_BER(txData,rxData);
        }
        //_BERprint(SNR_STEP, SNRdB, BER);

    }
    /*送受信データのプロット*/
    gp = popen("gnuplot -persist","w");//パイプを開き，gnuplotを立ち上げる
    //fprintf(gp, "set multiplot\n");//マルチプロットモード
    fprintf(gp,"set xrange [0:%d]\n",L);//範囲の指定
    fprintf(gp,"set yrange [0:1]\n");//範囲の指定
    fprintf(gp, "set xlabel \"data length\"\n"); // ラベル表示
    fprintf(gp, "set ylabel \"data\"\n");
    fprintf(gp,"plot '-' with lines linetype 1 title \"txData\",'-' with lines linetype 3 title \"rxData\" \n");
    /* 送信データ */
    for(i=0;i<L;++i){
        fprintf(gp,"%d\t%d\n",x[i],*(txData.dat+i));//送信データをプロット
    }
    fprintf(gp,"e\n");
    /* 受信データ */
    for(i=0;i<L;++i){
        fprintf(gp,"%d\t%d\n",x[i],*(rxData.dat+i));//受信データをプロット
    }
    fprintf(gp,"e\n");

    //fprintf(gp,"set nomultiplot\n");// マルチプロットモード終了
    fprintf(gp, "exit\n");//gnuplotの終了
    fflush(gp);
    pclose(gp);//パイプを閉じる


    /* 別ウインドウ */
    /* BERを表示 */
    gp2 = popen("gnuplot -persist","w");//パイプを開き，gnuplotを立ち上げる
    fprintf(gp2,"set xrange [0:%d]\n",SNR_STEP);//範囲の指定
    fprintf(gp2,"set yrange [0:0.1]\n");//範囲の指定
    fprintf(gp2, "set xlabel \"SNR [dB]\"\n"); // ラベル表示
    fprintf(gp2, "set ylabel \"BER\"\n");
    /* BERデータ */
    fprintf(gp,"plot '-' with lines linetype 1 title \"BER\"\n");
    for(i=0;i<SNR_STEP;++i){
        fprintf(gp,"%f\t%f\n",*(SNRdB+i),*(BER+i));//BERをプロット
    }
    fprintf(gp2,"e\n");

    fprintf(gp2,"set nomultiplot\n");// マルチプロットモード終了
    fprintf(gp2, "exit\n");//gnuplotの終了
    fflush(gp2);
    pclose(gp2);//パイプを閉じる

    free(txData.dat);free(rxData.dat);
    free(txSymbol.sig);free(rxSymbol.sig);
    free(noise.sig);
    return 0;
}