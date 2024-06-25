package sly.gwas.mem;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.TreeMap;
/**
 * 存储维护PLINK格式的GWAS数据的类；
 */
public class MemPLINK extends Mem{
    byte [] chrs=null;
    char [][] genos=null;
    long [] poss=null;
    TreeMap<Integer,Boolean> mSamples=null;
    public MemPLINK(String filenameTPED,String filenameTFAM){
        super();
        System.out.printf("read tfam file\n");
        if(!fillM(filenameTFAM)) {
            System.out.printf("tfam file format error\n");
        }
        else {
            System.out.printf("read tped file\n");
            if(!fillN(filenameTPED)){
                System.out.printf("tped file format error\n");
            }
            else{
                System.out.printf("#samples:%d\t#snps:%d\nloading data\n",m,n);
                if(!fillMem(filenameTPED)){
                    System.out.printf("fillMem error\n");
                }
                else{
                    System.out.printf("loaded\n");
                    ready=true;
                }
            }
        }
    }

    /**
     * 根据tfam文件的内容，获取数据中的样本信息；
     * @param filenameTFam tfam文件路径；
     * @return 是否成功；
     */
    private final boolean fillM(String filenameTFam) {
        mSamples=new TreeMap<Integer,Boolean>();
        File file = new File(filenameTFam);
        List<String> lines = null;
        try {
            lines=FileUtils.readLines(file, "UTF-8");
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        for(String line:lines) {
            String[] ss=line.split("\t");
            if(ss[ss.length-1].equals("1")) {
                mSamples.put(m,false);
                m0++;
                m++;
            }
            else if(ss[ss.length-1].equals("2")) {
                mSamples.put(m, true);
                m1++;
                m++;
            }
            else{
                m++;
                System.out.printf("error: unsupported phenotype code\n");
                return false;
            }
        }
        return true;
    }

    /**
     * 根据tped文件的内容，获取SNP数目；
     * @param filenameTPED tped文件路径；
     * @return 是否成功；
     */
    private final boolean fillN(String filenameTPED) {
        n=0;
        LineIterator it=null;
        try {
            it=FileUtils.lineIterator(new File(filenameTPED), "UTF-8");
            while (it.hasNext()) {
                it.nextLine();
                n++;
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        } finally {
            if(it!=null) {
                try {
                    it.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return true;
    }

    /**
     * 根据tped文件的内容，载入GWAS数据到内存；
     * @param filenameTPED tped文件路径；
     * @return 是否成功；
     */
    private final boolean fillMem(String filenameTPED) {
        LineIterator it=null;
        try {
            it=FileUtils.lineIterator(new File(filenameTPED), "UTF-8");
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        l1=(int)Math.ceil(((double)m1)/SIZE);
        l0=(int)Math.ceil(((double)m0)/SIZE);
        mem1=new long[n][3][l1];
        if(mem1==null)
            return false;
        mem0=new long[n][3][l0];
        if(mem0==null)
            return false;
        names=new String[n];
        if(names==null)
            return false;
        chrs=new byte[n];
        if(chrs==null)
            return false;
        poss=new long[n];
        if(poss==null)
            return false;
        genos=new char[n][2];
        if(genos==null)
            return false;
        int indexOfVector,indexInVector;
        char[] geno=new char[2];
        int[] cGeno=new int[2];
        int gg=0;
        char c=0;
        for(int i=0;i<n;i++){
            String[] ss=it.nextLine().split("\t");
            if(ss[0].equals("X")||ss[0].equals("x")) {
                chrs[i]=23;
            }
            else if(ss[0].equals("Y")||ss[0].equals("y")) {
                chrs[i]=24;
            }
            else {
                try {
                    chrs[i]=Byte.parseByte(ss[0]);
                }
                catch(java.lang.NumberFormatException e) {
                    chrs[i]=0;
                    System.out.printf("WARN: Can't get the chromosome number of %s(%s)\n",ss[0],ss[1]);
                }
            }
            names[i]=ss[1];
            poss[i]=Long.parseLong(ss[3]);
            cGeno[0]=0;
            cGeno[1]=0;
            geno[0]='X';
            geno[1]='X';
            for(int j=4;j<ss.length;j++) {
                c=ss[j].charAt(0);
                if(cGeno[0]==0) {
                    geno[0]=c;
                    cGeno[0]=1;
                }
                else {
                    if(c==geno[0]) {
                        cGeno[0]++;
                    }
                    else {
                        if(cGeno[1]==0) {
                            geno[1]=c;
                            cGeno[1]=1;
                        }
                        else {
                            if(c==geno[1]) {
                                cGeno[1]++;
                            }
                            else {
                                System.out.printf("genotype error\n");
                                return false;
                            }
                        }
                    }
                }
                c=ss[j].charAt(2);
                if(cGeno[0]==0) {
                    geno[0]=c;
                    cGeno[0]=1;
                }
                else {
                    if(c==geno[0]) {
                        cGeno[0]++;
                    }
                    else {
                        if(cGeno[1]==0) {
                            geno[1]=c;
                            cGeno[1]=1;
                        }
                        else {
                            if(c==geno[1]) {
                                cGeno[1]++;
                            }
                            else {
                                System.out.printf("genotype error\n");
                                return false;
                            }
                        }
                    }
                }
            }
            if(cGeno[0]<cGeno[1]){
                c=geno[0];
                geno[0]=geno[1];
                geno[1]=c;
            }
            genos[i][0]=geno[0];
            genos[i][1]=geno[1];
            for(int j=0, j0=0 ,j1=0;j<m;j++){
                c=ss[j+4].charAt(0);
                gg=0;
                if(c==geno[1]){
                    gg++;
                }
                c=ss[j+4].charAt(2);
                if(c==geno[1]){
                    gg++;
                }
                if(mSamples.get(j)) {
                    indexOfVector=j1/SIZE;
                    indexInVector=j1-SIZE*indexOfVector;
                    mem1[i][gg][indexOfVector]|=ONE>>>indexInVector;
                    j1++;
                }
                else {
                    indexOfVector=j0/SIZE;
                    indexInVector=j0-SIZE*indexOfVector;
                    mem0[i][gg][indexOfVector]|=ONE>>>indexInVector;
                    j0++;
                }
            }
        }
        return true;
    }
}
