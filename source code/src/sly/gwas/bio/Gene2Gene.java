package sly.gwas.bio;
import org.apache.commons.io.FileUtils;
import sly.gwas.mem.Mem;
import sly.gwas.mem.MemPLINK;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
/**
 * 只考虑基因与基因之间的直接关系；
 */
public class Gene2Gene {
    ArrayList<ArrayList<Integer>> gene2Gene=null;
    boolean ready=false;
    public Gene2Gene(String filenameGene2Gene, TreeMap<String,Integer> gene2Index){
        this.ready=false;
        this.gene2Gene=new ArrayList<>(gene2Index.size());
        for(int i=0;i<gene2Index.size();i++){
            this.gene2Gene.add(null);
        }
        if(!loadFromFile(filenameGene2Gene,gene2Index)){
            System.out.println("Gene2Gene::load file error");
        }
        else{
            ready=true;
        }
    }
    private boolean loadFromFile(String filenameGene2Gene, TreeMap<String,Integer> gene2Index){
        File file = new File(filenameGene2Gene);
        List<String> lines;
        try {
            lines = FileUtils.readLines(file, "UTF-8");
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        for(String s:lines){
            String[] ss=s.split("\t");
            int a=gene2Index.getOrDefault(ss[0],-1);
            int b=gene2Index.getOrDefault(ss[1],-1);
            if(a!=-1&&b!=-1){
                add2Gene2Gene(a,b);
                add2Gene2Gene(b,a);
            }
        }
        return true;
    }
    private void add2Gene2Gene(int a,int b){
        ArrayList<Integer> listA=this.gene2Gene.get(a);
        if(listA==null){
            listA=new ArrayList<>();
            listA.add(b);
            this.gene2Gene.set(a,listA);
        }
        else{
            if(!listA.contains(b))
                listA.add(b);
        }
    }
    /**
     * 测试函数；
     */
    public void test(ArrayList<String> genes){
        for(int i=0;i<gene2Gene.size();i++){
            System.out.printf("%s:",genes.get(i));
            ArrayList<Integer> a=gene2Gene.get(i);
            if(a!=null) {
                for (int j = 0; j < a.size();j++){
                    System.out.printf("\t\t%s",genes.get(a.get(j)));
                }
            }
            else{
                System.out.print("\tnull");
            }
            System.out.println();
        }
    }
    public static void main(String[] args){
        String filenameTFAM="E:\\workspace\\data\\gwas\\WTCCC20181025\\backup5\\bd_gwas.tfam";
        String filenameTPED="E:\\workspace\\data\\gwas\\WTCCC20181025\\backup5\\bd_gwas.tped";
        String filenameSNP2Gene="snps_after_2.txt";
        String filenameGene2Gene="gene_pairs_after_3.txt";
        Mem mem=new MemPLINK(filenameTPED,filenameTFAM);
        SNP2Gene s=new SNP2Gene(filenameSNP2Gene,mem);
        Gene2Gene g=new Gene2Gene(filenameGene2Gene,s.gene2Index);
        g.test(s.genes);
    }

    public boolean isReady() {
        return ready;
    }
    public final ArrayList<Integer> getGenesOfGene(int gene) {
        return this.gene2Gene.get(gene);
    }
}
