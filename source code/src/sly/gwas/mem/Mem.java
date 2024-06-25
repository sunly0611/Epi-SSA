package sly.gwas.mem;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import sly.gwas.Funcs;
import sly.gwas.Logger;
import sly.gwas.snp.SNPCom;
import sly.gwas.snp.SNPComWithG;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;

/**
 * 对GWAS数据的封装，GWAS数据的来源，目前有模拟数据和真实数据，对这两种数据格式做统一封装；
 */
public abstract class Mem{
    public String[] names=null;
    public int n;
    static int SIZE=Long.SIZE;
    static long ONE=0x8000000000000000l;
    long[][][] mem1=null;
    long[][][] mem0=null;
    int m,m1,m0,l1,l0;
    boolean ready=false;
    private static int[] BITS=_INITIAL_BITS();
    public Mem(){
        ready=false;
        BITS=new int[0x10000];
        for(int i=0;i<0x10000;i++){
            BITS[i]=bitCount(i);
        }
    }
    /**
     * 获取SNP组合的列联表；
     * @param snps SNP组合；
     * @return 列联表；
     */
    public final int[][] getTable(int[] snps){
        int lg=1;
        for(int i=0;i<snps.length;i++){
            lg=lg*3;
        }
        int [][] c=new int[2][lg];
        long[][] a=new long[lg][l0];
        for(int i=0;i<lg;){
            System.arraycopy(mem0[snps[0]][0],0,a[i++],0,l0);
            System.arraycopy(mem0[snps[0]][1],0,a[i++],0,l0);
            System.arraycopy(mem0[snps[0]][2],0,a[i++],0,l0);
        }
        for(int i=1,step=3;i<snps.length;i++,step*=3){
            for(int j=0;j<lg;j++){
                for(int k=0;k<l0;k++){
                    a[j][k]&=mem0[snps[i]][(j/step)%3][k];
                }
            }
        }
        for(int i=0;i<lg;i++){
            c[0][i]=popCount(a[i]);
        }
        a=new long[lg][l1];
        for(int i=0;i<lg;){
            System.arraycopy(mem1[snps[0]][0],0,a[i++],0,l1);
            System.arraycopy(mem1[snps[0]][1],0,a[i++],0,l1);
            System.arraycopy(mem1[snps[0]][2],0,a[i++],0,l1);
        }
        for(int i=1,step=3;i<snps.length;i++,step*=3){
            for(int j=0;j<lg;j++){
                for(int k=0;k<l1;k++){
                    a[j][k]&=mem1[snps[i]][(j/step)%3][k];
                }
            }
        }
        for(int i=0;i<lg;i++){
            c[1][i]=popCount(a[i]);
        }
        return c;
    }

    /**
     * 中间函数；
     * @param a
     * @return
     */
    private final static int popCount(long[] a) {
        int r=0;
        for(int i=0;i<a.length;i++){
            r=r+POP_COUNT(a[i]);
        }
        return r;
    }

    /**
     * 中间函数；
     * @param i
     * @return
     */
    private final static int POP_COUNT(long i ){
        return BITS[(int)(i&0xFFFF)] + BITS[(int)((i>>>16)&0xFFFF)] + BITS[(int)((i>>>32)&0xFFFF)] + BITS[(int)((i>>>48)&0xFFFF)];
    }

    /**
     * 中间函数；
     * @return
     */
    private final static int[] _INITIAL_BITS(){
        int[] r=new int[0x10000];
        for(int i=0;i<r.length;i++){
            r[i]=bitCount(i);
        }
        return r;
    }

    /**
     * 中间函数；
     * @param i
     * @return
     */
    private static final int bitCount(long i) {
        i = i - ((i >>> 1) & 0x5555555555555555l);
        i = (i & 0x3333333333333333l) + ((i >>> 2) & 0x3333333333333333l);
        i = (i + (i >>> 4)) & 0x0f0f0f0f0f0f0f0fl;
        i = i + (i >>> 8);
        i = i + (i >>> 16);
        i = i + (i >>> 32);
        return (int)i & 0x7f;
    }

    /**
     * 计算指定SNP组合的k2值；
     * @param indexes SNP组合；
     * @return k2值；
     */
    public final double computeK2(int [] indexes) {
        int [][] ca=getTable(indexes);
        int lg=ca[0].length;
        double k2=0;
        for(int i=0;i<lg;i++){
            double t=0;
            int ca0=ca[0][i];
            int ca1=ca[1][i];
            if(ca0==0&&ca1==0){
                continue;
            }
            else if(ca0>ca1){
                for(int j=2;j<=ca1;j++){
                    t+=Math.log(j);
                }
                for(int j=ca0+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            else{
                for(int j=2;j<=ca0;j++){
                    t+=Math.log(j);
                }
                for(int j=ca1+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            k2=k2+t;
        }
        k2=-k2;
        return k2;
    }
    /**
     * 计算指定SNP组合的k2值；
     * @param ca SNP组合的列联表；
     * @return k2值；
     */
    public final double computeK2(int [][] ca) {
        int lg=ca[0].length;
        double k2=0;
        for(int i=0;i<lg;i++){
            double t=0;
            int ca0=ca[0][i];
            int ca1=ca[1][i];
            if(ca0==0&&ca1==0){
                continue;
            }
            else if(ca0>ca1){
                for(int j=2;j<=ca1;j++){
                    t+=Math.log(j);
                }
                for(int j=ca0+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            else{
                for(int j=2;j<=ca0;j++){
                    t+=Math.log(j);
                }
                for(int j=ca1+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            k2=k2+t;
        }
        k2=-k2;
        return k2;
    }
    /**
     * 计算SNP组合的g值；
     * @param indexes SNP组合；
     * @return g值；
     */
    public final double computeG(int [] indexes) {
        int [][] ca=getTable(indexes);
        int lg=ca[0].length;
        double g=0;
        double p0=((double)m0)/m;
        double df=0;
        int c0,c1;
        double e;
        for(int i=0;i<lg;i++){
            c0=ca[0][i];
            c1=ca[1][i];
            if(c0!=0){
                e=(c0+c1)*p0;
                g+=c0*Math.log(c0/e);
            }
            if(c1!=0){
                e=(c0+c1)*(1-p0);
                g+=c1*Math.log(c1/e);
            }
            if(c0!=0||c1!=0)
                df++;
        }
        ChiSquaredDistribution d=new ChiSquaredDistribution(df-1);
        g=1-d.cumulativeProbability(g);
        if(Double.isNaN(g))
            g=1;
        return g;
    }

    /**
     * 计算SNP组合的g值；
     * @param ca SNP组合的列联表；
     * @return g值；
     */
    public final double computeG(int [][] ca) {
        int lg=ca[0].length;
        double g=0;
        double p0=((double)m0)/m;
        double df=0;
        int c0,c1;
        double e;
        for(int i=0;i<lg;i++){
            c0=ca[0][i];
            c1=ca[1][i];
            if(c0!=0){
                e=(c0+c1)*p0;
                g+=c0*Math.log(c0/e);
            }
            if(c1!=0){
                e=(c0+c1)*(1-p0);
                g+=c1*Math.log(c1/e);
            }
            if(c0!=0||c1!=0)
                df++;
        }
        ChiSquaredDistribution d=new ChiSquaredDistribution(df-1);
        g=1-d.cumulativeProbability(g);
        if(Double.isNaN(g))
            g=1;
        return g;
    }
    /**
     * 计算SNP组合的CE值；
     * @param indexes SNP组合；
     * @return ce；
     */
    public final double computeCe(int[] indexes){
        int [][] ca=getTable(indexes);
        return computeCe(ca);
    }

    /**
     * 计算SNP组合的CE值；
     * @param ca SNP组合的列联表；
     * @return ce；
     */
    public final double computeCe(int[][] ca){
        int lg=ca[0].length;
        double[][] d=new double[2][lg];
        double[] px=new double[lg];
        for(int i=0;i<2;i++){
            for(int j=0;j<lg;j++){
                d[i][j]=((double)ca[i][j])/m;
            }
        }
        for(int i=0;i<lg;i++){
            px[i]=d[0][i]+d[1][i];
        }
        double ce=0;
        for(int i=0;i<2;i++){
            for(int j=0;j<lg;j++){
                if(d[i][j]!=0)
                    ce+=d[i][j]*Math.log(d[i][j]/px[j]);
            }
        }
        ce=-ce;
        return ce;
    }

    /**
     * 计算SNP组合的AdjustCE值；
     * @param indexes SNP组合；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return 修正后的ce；
     */
    public final double computeAdjustCe(int[] indexes,int maxL){
        int [][] ca=getTable(indexes);
        return computeAdjustCe(ca,maxL);
    }

    /**
     * 计算SNP组合的AdjustCE值；
     * @param ca SNP组合的列联表；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return 修正后的ce；
     */
    public final double computeAdjustCe(int[][] ca,int maxL){
        int lg=ca[0].length;
        double ce=0;
        double p0=0;
        double p1=0;
        double p=0;
        double v=0;
        ArrayList<ValueIndex> vis=new ArrayList<>();
        for(int i=0;i<lg;i++){
            if(ca[0][i]!=0||ca[1][i]!=0) {
                p0 = ((double) ca[0][i]) / m;
                p1 = ((double) ca[1][i]) / m;
                p = p0 + p1;
                v = 0;
                if(p0!=0)
                    v=v-p0*Math.log(p0/p);
                if(p1!=0)
                    v=v-p1*Math.log(p1/p);
                vis.add(new ValueIndex(v,i));
            }
        }
        if(vis.size()>maxL){
            vis.sort(ValueIndexComparator);
            for(int i=0;i<maxL-1;i++){
                ce+=vis.get(i).v;
            }
            int ca0=0;
            int ca1=0;
            for(int i=maxL-1;i<vis.size();i++){
                ValueIndex vi=vis.get(i);
                ca0+=ca[0][vi.i];
                ca1+=ca[1][vi.i];
            }
            p0 = ((double) ca0) / m;
            p1 = ((double) ca1) / m;
            p = p0 + p1;
            if(p0!=0)
                ce=ce-p0*Math.log(p0/p);
            if(p1!=0)
                ce=ce-p1*Math.log(p1/p);
        }
        else{
            for(ValueIndex vi:vis){
                ce+=vi.v;
            }
        }

        return ce;
    }

    /**
     * 计算gini的值；
     * @param ca SNP组合的列联表；
     * @return gini值；
     */
    public final double computeGini(int [][] ca){
        int lg=ca[0].length;
        double gini=0;
        for(int i=0;i<lg;i++){
            double a=(double)(ca[0][i]+ca[1][i]);
            if(a!=0){
                double p=ca[0][i]/a;
                gini+=a/m*2*p*(1-p);
            }
        }
        return gini;
    }

    /**
     * 计算gini的值；
     * @param indexes SNP组合；
     * @return gini值；
     */
    public final double computeGini(int[] indexes){
        int [][] ca=getTable(indexes);
        return computeGini(ca);
    }
    /**
     * 计算AdjustGini的值；
     * 1.为了消除对列联表长度的偏性，在通过累加计算gini值时，只考虑最小的maxL-1个，对于剩余的部分，对列联表进行样本上的合并，合并为第maxL个单元，计算gini值，实际上，就是对列联表的合并；
     * 2.对于修正后的列联表，如常计算gini的值；
     * @param ca SNP组合的列联表；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return AdjustGini的值；
     */
    public final double computeAdjustGini(int[][] ca,int maxL){
        int lg=ca[0].length;
        ArrayList<ValueIndex> vis=new ArrayList<>();
        double gini=0;
        for(int i=0;i<lg;i++){
            double a=(double)(ca[0][i]+ca[1][i]);
            if(a!=0){
                double p=ca[0][i]/a;
                vis.add(new ValueIndex(a/m*2*p*(1-p),i));
            }
        }
        if(vis.size()>maxL){
            vis.sort(ValueIndexComparator);
            for(int i=0;i<maxL-1;i++){
                gini+=vis.get(i).v;
            }
            int ca0=0;
            int ca1=0;
            for(int i=maxL-1;i<vis.size();i++){
                ValueIndex vi=vis.get(i);
                ca0+=ca[0][vi.i];
                ca1+=ca[1][vi.i];
            }
            double a=(double)(ca0+ca1);
            double p=ca0/a;
            gini+=a/m*2*p*(1-p);
        }
        else{
            for(ValueIndex vi:vis){
                gini+=vi.v;
            }
        }
        return gini;
    }

    /**
     * 计算SNP组合的AdjustGini的值；
     * 根据评测经验，在评测一个上位性组合的gini值时，对于一个与疾病无关的上位性组合，如果将样本分的过细，或者说，列联表过宽，则即使一个与疾病统计上无关的上位性组合也会产生较强的gini值；
     * 例如：
     * gini(sc)的值与基于sc构建的列联表的长度之间有很大的偏性。
     * @param indexes 一个snp组合；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return 修正后的gini值；
     */
    public final double computeAdjustGini(int[] indexes,int maxL){
        int[][] ca=getTable(indexes);
        return computeAdjustGini(ca,maxL);
    }

    /**
     * 计算LR的值；
     * 来自于文献；
     * 原文中LR值越大，关联性越强，因此我在这里实际计算的是NLR，取了负值；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param ca 列联表；
     * @return LR值；
     */
    public final double computeLR(int[][] ca){
        int lg=ca[0].length;
        double lr=0;
        double p0=((double)m0)/m;
        int c0,c1;
        double e;
        for(int i=0;i<lg;i++){
            c0=ca[0][i];
            c1=ca[1][i];
            if(c0!=0){
                e=(c0+c1)*p0;
                lr+=c0*Math.log(c0/e);
            }
            if(c1!=0){
                e=(c0+c1)*(1-p0);
                lr+=c1*Math.log(c1/e);
            }
        }
        return -2*lr;
    }
    /**
     * 计算LR的值；来自于文献；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param indexes snp组合；
     * @return LR值；
     */
    public final double computeLR(int[] indexes){
        int [][] ca=getTable(indexes);
        return computeLR(ca);
    }

    /**
     * 计算JS的值；
     * 来自于文献；
     * 原文中JS值越大，关联性越强，因此我在这里实际计算的是NJS，取了负值；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param ca 列联表；
     * @return JS值；
     */
    public final double computeJS(int[][] ca){
        int lg=ca[0].length;
        double js=0;
        for(int i=0;i<lg;i++){
            int c=ca[0][i]+ca[1][i];
            if(ca[0][i]>0) {
                double t = ca[0][i];
                t = t / c;
                js+=Math.log(2*t)*t;
            }
            if(ca[1][i]>0) {
                double t = ca[1][i];
                t = t / c;
                js+=Math.log(2*t)*t;
            }
        }
        return -js/2;
    }

    /**
     * 计算JS的值；来自于文献；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param indexes snp组合；
     * @return LR值；
     */
    public final double computeJS(int[] indexes){
        int [][] ca=getTable(indexes);
        return computeJS(ca);
    }

    //计算NJ的前置；
    private double [] d=null;
    private final double computeD(int[] indexes){
        double d=0;
        int [][] ca=getTable(indexes);
        int lg=ca[0].length;
        for(int i=0;i<lg;i++){
            d+=Math.abs(ca[0][i]-ca[1][i]);
        }
        return d;
    }
    /**
     * 计算NJ的值；来自于文献；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param ca 列联表；
     * @return LR值；
     */
    public final double computeNJ(int[][] ca,int[] indexes){
        if(d==null){
            d=new double[this.n];
            for(int i=0;i<this.n;i++){
                d[i]=computeD(new int[]{i});
            }
        }
        double D=0;
        int lg=ca[0].length;
        for(int i=0;i<lg;i++){
            D+=Math.abs(ca[0][i]-ca[1][i]);
        }
        double ND=0;
        for(int i:indexes){
            ND+=d[i];
        }
        ND=ND/D;
        double JE=0;
        for(int i=0;i<lg;i++){
            if(ca[0][i]>0) {
                double t=ca[0][i];
                t=t/m0;
                JE += Math.log(t) * t;
            }
        }
        JE=-JE;
        return ND/JE;
    }
    /**
     * 计算NJ的值；来自于文献；
     * Multi-population harmony search algorithm for the
     * detection of high-order SNP interactions
     * @param indexes snp组合；
     * @return LR值；
     */
    public final double computeNJ(int[] indexes){
        int[][] ca=getTable(indexes);
        return computeNJ(ca,indexes);
    }
    public final static Comparator<ValueIndex> ValueIndexComparator=new Comparator<ValueIndex>() {
        //从小到大排序
        @Override
        public int compare(ValueIndex o1, ValueIndex o2) {
            if(o1.v<o2.v){
                return 1;
            }
            else if(o1.v>o2.v){
                return -1;
            }
            else{
                return 0;
            }
        }
    };

    class ValueIndex{
        double v=0;
        int i=0;
        public ValueIndex(double v,int i){
            this.v=v;
            this.i=i;
        }
    }
    /**
     * 计算k2值的函数，尝试修正由于分裂太细而导致的对k2值的偏性；
     * 将列联表修正到固定的长度，以防评价指标对细分列联表的偏性；
     * 通过k2的值进行修复；
     * @param ca 原始列联表；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return 修正的k2值；
     */
    public final double computeAdjustK2(int[][] ca,int maxL){
        int lg=ca[0].length;
        double k2=0;
        ArrayList<ValueIndex> ks=new ArrayList<>();
        for(int i=0;i<lg;i++){
            double t=0;
            int ca0=ca[0][i];
            int ca1=ca[1][i];
            if(ca0==0&&ca1==0){
                continue;
            }
            else if(ca0>ca1){
                for(int j=2;j<=ca1;j++){
                    t+=Math.log(j);
                }
                for(int j=ca0+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            else{
                for(int j=2;j<=ca0;j++){
                    t+=Math.log(j);
                }
                for(int j=ca1+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            ks.add(new ValueIndex(-t,i));
        }
        if(ks.size()>maxL){
            ks.sort(ValueIndexComparator);
            for(int i=0;i<maxL-1;i++){
                k2+=ks.get(i).v;
            }
            int ca0=0;
            int ca1=0;
            for(int i=maxL-1;i<ks.size();i++){
                ValueIndex vi=ks.get(i);
                ca0+=ca[0][vi.i];
                ca1+=ca[1][vi.i];
            }
            double t=0;
            if(ca0>ca1){
                for(int j=2;j<=ca1;j++){
                    t+=Math.log(j);
                }
                for(int j=ca0+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            else{
                for(int j=2;j<=ca0;j++){
                    t+=Math.log(j);
                }
                for(int j=ca1+1;j<=ca0+ca1+1;j++){
                    t-=Math.log(j);
                }
            }
            k2=k2-t;
        }
        else{
            for(ValueIndex vi:ks){
                k2+=vi.v;
            }
        }
        return k2;

    }

    /**
     *
     * @param indexes snp组合；
     * @param maxL 合并列联表时，合并后的列联表的最大长度；
     * @return 修正的k2值；
     */
    public final double computeAdjustK2(int[] indexes,int maxL){
        int [][] ca=getTable(indexes);
        return computeAdjustK2(ca,maxL);
    }

    /**
     * 推荐order；
     * @return 推荐的order；
     */
    public final int recommendOrder() {
        double t=Math.log(Math.min(m0,m1))-0.5;
        return (int)Math.floor(t);
        /*
        double t=Math.log(Math.min(m0,m1)/5)/Math.log(3);
        return (int)Math.floor(t);
        */
    }

    public final int recommendMaxL(){
        return Math.min(m0,m1)/10;
    }

    /**
     * 测试函数；
     */
    public final void testMem() {
        for(int g=0;g<3;g++) {
            for(int i=0;i<n;i++) {
                System.out.printf("mem0 g%d snp%d: ",g,i);
                for(int j=0;j<l0;j++) {
                    System.out.printf("%s",Long.toBinaryString(mem0[i][g][j]));
                }
                System.out.println();
            }
        }
        for(int g=0;g<3;g++) {
            for(int i=0;i<n;i++) {
                System.out.printf("mem1 g%d snp%d: ",g,i);
                for(int j=0;j<l1;j++) {
                    System.out.printf("%s",Long.toBinaryString(mem1[i][g][j]));
                }
                System.out.println();
            }
        }
    }

    /**
     * 测试函数；
     */
    public final void testTable() {
        for(int i=0;i<n;i++) {
            int [] x= {i};
            int[][] ca=getTable(x);
            int c0=0;
            int c1=0;
            for(int j=0;j<ca[0].length;j++) {
                c0+=ca[0][j];
                c1+=ca[1][j];
            }
            System.out.printf("%s: %d %d %d\n",names[i],c0,c1,c0+c1);
        }
    }

    /**
     * 获取一阶g的值；
     */
    public final double[] computeG1() {
        double[] r=new double[n];
        int[] x=new int[1];
        for(int i=0;i<n;i++){
            x[0]=i;
            r[i]=computeG(x);
        }
        return r;
    }

    /**
     * 获取SNP名称与内存下标的映射；
     * @return 映射；
     */
    public final TreeMap<String,Integer> reversedNames(){
        TreeMap<String,Integer> r=new TreeMap<String,Integer>();
        for(int i=0;i<n;i++) {
            r.put(names[i], i);
        }
        return r;
    }

    /**
     * 测试函数；
     * @param x snp组合；
     */
    public void testTable(int[] x){
            int[][] ca=getTable(x);
        System.out.println(Arrays.toString(ca[0]));
        System.out.println(Arrays.toString(ca[1]));
        int zero=0;
        for(int i=0;i<ca[0].length;i++){
            if(ca[0][i]==0&&ca[1][i]==0){
                zero++;
            }
        }
        System.out.println(String.format("not zero: %d",ca[0].length-zero));
        System.out.println(String.format("zero: %d",zero));
        System.out.println(String.format("k2: %f",computeK2(ca)));
        System.out.println(String.format("g: %f",computeG(ca)));
        System.out.println(String.format("ce: %f",computeCe(ca)));
        System.out.println(String.format("gini: %f",computeGini(ca)));
    }

    /**
     * 在x上检测上位性，基于k2指标；
     * @param x 输入的SNP组合；
     * @return 在SNP组合上检测到的上位性组合，如果不存在，返回null。
     */
    public final int[] detectEpisBasedOnK2(int[] x){
        //Funcs.testSolution(x);
        Logger.debug(String.format("detect on %s\n",Arrays.toString(x)));
        boolean successReduced=true;
        double k2=computeK2(x);
        while(successReduced&&x.length>=2) {
            successReduced=false;
            int[] xReduced=new int[x.length-1];
            double k2R=0;
            for (int i = 0; i < x.length; i++) {
                for (int j = 0; j < i; j++) {
                    xReduced[j] = x[j];
                }
                for (int j = i; j < xReduced.length; j++) {
                    xReduced[j] = x[j + 1];
                }
                k2R = this.computeK2(xReduced);
                if (k2R <= k2) {
                    successReduced = true;
                    x = xReduced;
                    k2 = k2R;
                    break;
                }
            }
        }
        if(x.length>=2) {
            return x;
        }
        else{
            return null;
        }
    }

    /**
     * 在x上检测上位性，基于k2指标；
     * 与之前的策略不同之处在于，不再那么信任k2指标了，每次缩减，都是找衰减最大的SNP进行缩减；
     * @param x 输入的SNP组合；
     * @return 在SNP组合上检测到的上位性组合，如果不存在，返回null。
     */
    public final int[] detectEpisBasedOnMaxReducedK2(int[] x){
        //Funcs.testSolution(x);
        Logger.debug(String.format("detect on %s\n",Arrays.toString(x)));
        double k2=computeK2(x);
        while(x.length>=2) {
            int indexBest=-1;
            double k2RBest=Double.MAX_VALUE;
            int[] xReduced=new int[x.length-1];
            for (int i = 0; i < x.length; i++) {
                for (int j = 0; j < i; j++) {
                    xReduced[j] = x[j];
                }
                for (int j = i; j < xReduced.length; j++) {
                    xReduced[j] = x[j + 1];
                }
                double k2R = this.computeK2(xReduced);
                if (k2R <= k2 && k2R<k2RBest) {
                    indexBest=i;
                    k2RBest=k2R;
                }
            }
            if(indexBest==-1) {
                break;
            }
            else{
                for (int j = 0; j < indexBest; j++) {
                    xReduced[j] = x[j];
                }
                for (int j = indexBest; j < xReduced.length; j++) {
                    xReduced[j] = x[j + 1];
                }
                x=xReduced;
                k2 = k2RBest;
            }

        }
        if(x.length>=2) {
            return x;
        }
        else{
            return null;
        }
    }

    /**
     * 在x上检测上位性，基于adjust k2指标；
     * @param x 输入的SNP组合；
     * @return 在SNP组合上检测到的上位性组合，如果不存在，返回null。
     */
    public final int[] detectEpisBasedOnAdjustK2(int[] x,int maxL){
        //Funcs.testSolution(x);
        Logger.debug(String.format("detect on %s\n",Arrays.toString(x)));
        boolean successReduced=true;
        double k2=computeAdjustK2(x,maxL);
        while(successReduced&&x.length>=2) {
            successReduced=false;
            int[] xReduced=new int[x.length-1];
            double k2R=0;
            for (int i = 0; i < x.length; i++) {
                for (int j = 0; j < i; j++) {
                    xReduced[j] = x[j];
                }
                for (int j = i; j < xReduced.length; j++) {
                    xReduced[j] = x[j + 1];
                }
                k2R = this.computeAdjustK2(xReduced,maxL);
                if (k2R <= k2) {
                    successReduced = true;
                    x = xReduced;
                    k2 = k2R;
                    break;
                }
            }
        }
        if(x.length>=2) {
            return x;
        }
        else{
            return null;
        }
    }

    /**
     *
     * @return
     */
    public final boolean isReady(){
        return ready;
    }

    public String getName(int i){
        return this.names[i];
    }
    final public static void main(String[] args) {
        String pathIn="E:\\workspace\\data\\gwas\\simulated_data\\DNME01_1600_1000\\7_750.txt";
        Mem mem=new MemSimulated(pathIn);
        int[] a=new int[]{7,124,454,750};
        int[] b=new int[] {75,563,824,842};
        mem.testTable(a);
        mem.testTable(b);
    }
}
