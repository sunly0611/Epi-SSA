package sly.gwas;

import java.util.Arrays;
import java.util.Random;

public final class Funcs {
    /**
     * 向x中添加元素，升序；
     * @param x 当前的数组，本身也是升序的，长度为pos-1；
     * @param pos 添加value后，x数组的长度；
     * @param value 添加的值；
     * @return void。
     */
    public final static void insertAsc(int[] x, int pos,int value){
        for(int i=0;i<pos;i++){
            if(x[i]>value){
                for(int j=pos;j>i;j--){
                    x[j]=x[j-1];
                }
                x[i]=value;
                return;
            }
        }
        x[pos]=value;
    }


    /**
     * 从0到n-1的n个数里面，随机抽取k个，组成升序的数组，返回；
     * @param random 随机数生成对象；
     * @param n 随机数的上界；
     * @param k 返回的随机数数目；
     * @return 长度为k的数组，升序。
     */
    public final static int[] sample(Random random, int n, int k){
        int[] r=new int[k];
        for(int i=0;i<k;i++){
            insertAsc(random,r,n,i);
        }
        return r;
    }

    /**
     * 从0到n-1之间随机选择一个整数，插入到x的第i个位置上，不可以重复，且保持x有序；
     * @param random 随机数生成类对象；
     * @param x 待插入随机元素的数组；函数执行完毕后，插入一个新的随机元素；
     * @param n n；
     * @param i 插入元素之前，数组的长度，也是插入元素之后，数组中最后一个有元素的位置；
     * @return null；
     */
    public final static void insertAsc(Random random,int[] x, int n,int i){
        //Logger.debug(String.format("insert random in %d to %s\n",i,Arrays.toString(x)));
        int ti=random.nextInt(n-i);
        for(int j=0;j<i;j++){
            if(ti>=x[j]){
                ti=ti+1;
            }
            else{
                break;
            }
        }
        insertAsc(x,i,ti);
    }

    /**
     * SNP组合的加法运算，模拟a与b的加法运算，也是a向b移动的过程；
     * @param random 随机数生成类对象；
     * @param a 一个SNP组合；
     * @param b 一个SNP组合；
     * @return a+b的结果；
     */
    public static final int[] add(Random random,int[] a,int[] b){
        int n=a.length;
        int[] r=new int[n];
        //记录a中被抽取的元素的地址；
        boolean[] selectedInA=new boolean[n];
        Arrays.fill(selectedInA,false);
        int countInA=0;
        //记录b中被抽取的元素的地址；
        boolean[] selectedInB=new boolean[n];
        Arrays.fill(selectedInB,false);
        int countInB=0;
        //从a中选择一半的元素
        for(int i=0;i<n/2;i++){
            int t = random.nextInt(n - countInA);
            while (selectedInA[t]) {
                t++;
            }
            int selectedValue = a[t];
            countInA++;
            selectedInA[t] = true;
            for (int j = 0; j < n; j++) {
                if (b[j] == selectedValue) {
                    selectedInB[j] = true;
                    countInB++;
                    break;
                }
            }
            insertAsc(r, i, selectedValue);
        }
        //从b中选择一半的元素
        for(int i=n/2;i<n;i++){
            int t = random.nextInt(n - countInB);
            while (selectedInB[t]) {
                t++;
            }
            int selectedValue = b[t];
            countInB++;
            selectedInB[t] = true;
            insertAsc(r, i, selectedValue);
        }
        return r;
    }
    /**
     * 从0到n-1的n个数里面，以概率pss，随机抽取k个，组成升序的数组，返回；
     * @param random 随机数生成对象；
     * @param n 随机数的上界；
     * @param k 返回的随机数数目；
     * @param pss 0到n-1中每个元素被选中的概率；
     * @return 长度为k的数组，升序。
     */
    public final static int[] sample(Random random, int n, int k, double[] pss){
        int[] r=new int[k];
        for(int i=0;i<k;i++){
            insertAsc(random,r,n,i,pss);
        }
        return r;
    }

    /**
     * 从0到n-1之间，以概率随机选择一个整数，插入到x的第i个位置上，不可以重复，且保持x有序；
     * @param random 随机数生成类对象；
     * @param x 待插入随机元素的数组；函数执行完毕后，插入一个新的随机元素；
     * @param n n；
     * @param i 插入元素之前，数组的长度，也是插入元素之后，数组中最后一个有元素的位置；
     * @param pss 0到n-1中每个元素被选中的概率；
     * @return null；
     */
    private final static void insertAsc(Random random,int[] x, int n,int i,double[] pss){
        double sum=0;
        int ti = -1;
        //0->x[0]-1   x[j-1]+1->x[j]-1  x[i-1]+1->n-1
        if(i>0) {
            for (int j = 0; j < x[0]; j++) {
                sum += pss[j];
            }
            for (int j = 1; j < i; j++) {
                for (int k = x[j - 1] + 1; k < x[j]; j++) {
                    sum += pss[k];
                }
            }
            for (int j = x[i - 1] + 1; j < n; j++) {
                sum += pss[j];
            }
            double rSum = random.nextDouble() * sum;
            for (int j = 0; j < x[0]; j++) {
                rSum -= pss[j];
                if (rSum <= 0) {
                    ti = j;
                    break;
                }
            }
            if (ti == -1) {
                for (int j = 1; j < i; j++) {
                    for (int k = x[j - 1] + 1; k < x[j]; j++) {
                        rSum -= pss[k];
                        if (rSum <= 0) {
                            ti = k;
                            break;
                        }
                    }
                }
            }
            if (ti == -1) {
                for (int j = x[i - 1] + 1; j < n; j++) {
                    rSum -= pss[j];
                    if (rSum <= 0) {
                        ti = j;
                        break;
                    }
                }
            }
        }
        else{
            for (int j = 0; j < n; j++) {
                sum += pss[j];
            }
            double rSum = random.nextDouble() * sum;
            for (int j = 0; j < n; j++) {
                rSum -= pss[j];
                if (rSum <= 0) {
                    ti = j;
                    break;
                }
            }
        }
        if(ti==-1){
            Logger.error("Funcs.insertAsc ti==-1\n");
        }

        insertAsc(x,i,ti);

    }

    /**
     * 测试用，检测一个x数组是否内部增序；
     * @param x 待检查的数组；
     * @return 有错误返回false，否则，满足要求，返回true。
     */
    public final static boolean check(int [] x){
        for(int i=0;i<x.length-1;i++){
            if(x[i]>=x[i+1]){
                Logger.error("check error\n");
                return false;
            }
        }
        return true;
    }

    /**
     * 测试x中是否包含对应的SNP集合；
     * @param x 待检测的SNP组合；
     */
    public final static void testSolution(int[] x) {
        boolean i1=false;
        boolean i2=false;
        for(int i=0;i<x.length;i++){
            if(x[i]==146){
                i1=true;
            }
            if(x[i]==231){
                i2=true;
            }
        }
        if(i1&&i2){
            Logger.debug("testSolution in\n");
        }
    }
}
