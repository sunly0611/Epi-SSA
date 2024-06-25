package sly.gwas.bio;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import org.apache.commons.io.FileUtils;
import sly.gwas.mem.Mem;
import sly.gwas.mem.MemPLINK;
/**
 * 存储，维护SNP与基因之间关系的类；
 */
final public class SNP2Gene {
    int nGenes=-1;
    boolean ready=false;
    ArrayList<String> genes=null;
    TreeMap<String,Integer> gene2Index=null;
    ArrayList<ArrayList<Integer>> snp2Genes=null;
    ArrayList<ArrayList<Integer>> gene2Snps=null;
    public SNP2Gene(String filenameSNP2Gene, Mem mem){
        ready=false;
        this.snp2Genes=new ArrayList<ArrayList<Integer>>(mem.n);
        this.genes=new ArrayList<String>();
        this.gene2Index=new TreeMap<String,Integer>();
        for(int i=0;i<mem.n;i++){
            snp2Genes.add(null);
        }
        if(!loadSNP2Genes(filenameSNP2Gene,mem)){
            System.out.printf("SNP2Gene::loadSNP2Genes error\n");
            return;
        }
        ready=true;
    }
    final boolean loadSNP2Genes(String filenameSNP2Gene,Mem mem){
        TreeMap<String,Integer> maps=mem.reversedNames();
        File file = new File(filenameSNP2Gene);
        List<String> lines;
        try {
            lines = FileUtils.readLines(file, "UTF-8");
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        this.gene2Index.put("unknown",0);
        int iGene=1;
        genes=new ArrayList<String>();
        genes.add("unknown");
        for(String line:lines) {
            String[] ss=line.split("\t");
            int key=maps.getOrDefault(ss[0],-1);
            if(key!=-1) {
                ArrayList<Integer> a=new ArrayList<Integer>();
                for(int i=1;i<ss.length;i++) {
                    if(this.gene2Index.containsKey(ss[i])) {
                        a.add(this.gene2Index.get(ss[i]));
                    }
                    else {
                        this.gene2Index.put(ss[i],iGene);
                        a.add(iGene);
                        genes.add(ss[i]);
                        iGene++;
                    }
                }
                this.snp2Genes.set(key,a);
            }
        }
        for(int i=0;i<mem.n;i++) {
            if(this.snp2Genes.get(i)==null) {
                ArrayList<Integer> a=new ArrayList<Integer>();
                a.add(0);
                this.snp2Genes.set(i,a);
            }
        }
        this.nGenes=gene2Index.size();
        this.gene2Snps=new ArrayList<ArrayList<Integer>>(this.nGenes);
        for(int i=0;i<nGenes;i++) {
            gene2Snps.add(null);
        }
        Iterator<ArrayList<Integer>> it=snp2Genes.iterator();
        for(int i=0;it.hasNext();i++) {
            ArrayList<Integer> a=it.next();
            if(a!=null) {
                Iterator<Integer> iit=a.iterator();
                while(iit.hasNext()) {
                    int snp=iit.next();
                    ArrayList<Integer> l=gene2Snps.get(snp);
                    if(l==null) {
                        l=new ArrayList<Integer>();
                        l.add(i);
                        gene2Snps.set(snp,l);
                    }
                    else {
                        l.add(i);
                    }
                }
            }
        }
        return true;
    }
    final void test(Mem mem) {
        if(snp2Genes!=null) {
            System.out.println("snp 2 gene:");
            for(int i=0;i<mem.n;i++) {
                if(snp2Genes.get(i)!=null) {
                    System.out.print(mem.names[i]);
                    for(int gene:snp2Genes.get(i)) {
                        System.out.printf(",%s", genes.get(gene));
                    }
                    System.out.println();
                }
            }
        }
        if(gene2Snps!=null) {
            System.out.println("gene 2 snp:");
            for(int i=0;i<nGenes;i++) {
                if(gene2Snps.get(i)!=null) {
                    System.out.print(genes.get(i));
                    for(int snp:gene2Snps.get(i)) {
                        System.out.printf(",%s", mem.names[snp]);
                    }
                    System.out.println();
                }
            }
        }
    }
    final public static void main(String[] args) {
        String filenameTFAM="E:\\workspace\\data\\gwas\\WTCCC20181025\\backup5\\bd_gwas.tfam";
        String filenameTPED="E:\\workspace\\data\\gwas\\WTCCC20181025\\backup5\\bd_gwas.tped";
        String filenameSNP2Gene="snps_after_2.txt";
        Mem mem=new MemPLINK(filenameTPED,filenameTFAM);
        SNP2Gene s=new SNP2Gene(filenameSNP2Gene,mem);
        s.test(mem);
    }

    public final boolean isReady() {
        return ready;
    }

    public final TreeMap<String, Integer> getGene2Index() {
        return gene2Index;
    }
    public final ArrayList<Integer> getGenesOfSNP(int i){
        return snp2Genes.get(i);
    }
    public final ArrayList<Integer> getSNPsOfGene(int i){
        return gene2Snps.get(i);
    }
}
