package sly.gwas.test;

import sly.gwas.mem.Mem;
import sly.gwas.mem.MemSimulated;
import sly.gwas.snp.SNPComWithG;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

public abstract class Test {
    /**
     * 向x中添加元素，升序；
     * @param x 当前的数组，本身也是升序的，长度为pos-1；
     * @param pos 添加value后，x数组的长度；
     * @param value 添加的值；
     */
    private void insertAsc(int[] x, int pos,int value){
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
     * 随机生成一个SNP组合；
     * @param o 阶数；
     * @return 随机组合；
     */
    private int[] genRandomSNPCombination(Random random,Mem mem,int o){
        int[] r=new int[o];
        int ti=-1;
        for(int i=0;i<o;i++){
            ti=random.nextInt(mem.n-i);
            for(int j=0;j<i;j++){
                if(ti>=r[j]){
                    ti=ti+1;
                }
                else{
                    break;
                }
            }
            insertAsc(r,i,ti);
        }
        return r;
    }

    /**
     * 测试o阶SNP组合中，k2指标可以起到的作用，该函数不会进行全局搜索，而是随机测试，是针对阶数太高，测试过于耗时设计的；
     * @param o 阶数；
     * @param n 测试次数；
     * @param pathO 用于存储测试结果的文件路径；
     */
    public final void run(int o,int n,String pathIn,String pathO,int or,int seed) throws IOException {
        Random random=new Random(seed);
        ArrayList<SNPComWithG> snps=new ArrayList<>();
        Mem mem = new MemSimulated(pathIn);
        for(int i=0;i<n;i++) {
            int[] x = genRandomSNPCombination(random, mem, o);
            snps.add(new SNPComWithG(x, computeScore(mem,x)));
        }
        Collections.sort(snps, new Comparator<SNPComWithG>() {
            @Override
            public int compare(SNPComWithG o1, SNPComWithG o2) {
                return -o1.compareTo(o2);
            }
        });
        BufferedOutputStream bos=new BufferedOutputStream(new FileOutputStream(pathO));
        Iterator<SNPComWithG> it=snps.iterator();
        while(it.hasNext()){
            bos.write(it.next().getOStringForTest(mem,or).getBytes());
        }
        bos.flush();
        bos.close();
    }

    public abstract double computeScore(Mem mem,int[] index);

    public static final void main(String[] args) throws IOException {
        //String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME01_1600_100\\2_19.txt";
        //String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME01_1600_1000\\7_750.txt";
        //String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME05_1600_1000\\1_893.txt";
        //String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME12_1600_1000\\10_525.txt";
        //String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME06_1600_1000\\0_47.txt";
        String pathIn="E:\\workspace\\data\\gwas\\sim_data_20220613\\DME01_1600_1000\\168_950.txt";
        if(!new File("test").exists()){
            new File("test").mkdir();
        }
        int n=2000000;
        int seed=0;
        for(int o:new int[]{4,5,6}){
            System.out.println("test on "+o);
            int maxL=800/10;
            System.out.println("maxL "+maxL);

            new Test(){
                @Override
                public double computeScore(Mem mem, int[] index) {
                    return mem.computeAdjustGini(index,maxL);
                }
            }.run(o,n,pathIn,String.format("test/%d_ag.txt",o),2,seed);
            new Test(){
                @Override
                public double computeScore(Mem mem, int[] index) {
                    return mem.computeAdjustK2(index,maxL);
                }
            }.run(o,n,pathIn,String.format("test/%d_ak2.txt",o),2,seed);

            new Test(){
                @Override
                public double computeScore(Mem mem, int[] index) {
                    return mem.computeAdjustCe(index,maxL);
                }
            }.run(o,n,pathIn,String.format("test/%d_ace.txt",o),2,seed);
        }
    }
}
