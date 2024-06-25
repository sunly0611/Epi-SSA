package sly.gwas.algorithm;
import org.apache.commons.cli.*;
import sly.gwas.bio.Gene2Gene;
import sly.gwas.bio.SNP2Gene;
import sly.gwas.mem.Mem;
import sly.gwas.mem.MemPLINK;
import sly.gwas.mem.MemSimulated;
import sly.gwas.snp.SNPComWithG;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
/**
 * 用于维护程序参数；
 */
public class SHEIB {
    String pathIn="data.txt";
    String pathO="result.txt";
    String pathS2G=null;
    String pathAG=null;
    int o=-1;
    int maxGen=-1;
    double pb=0.8;
    double cG=0.05;
    double cGc=0.05;
    int rn=-1;
    int type=0;
    int seed=0;
    Mem mem=null;
    SNP2Gene snp2Gene=null;
    Gene2Gene gene2Gene=null;
    Random random=null;
    boolean ready=false;
    TreeSet<SNPComWithG> results=null;
    public SHEIB(String [] args){
        ready=false;
        Options options = new Options();
        options.addOption("pathIn", true, "path of the GWAS data, default 'data.txt'.");
        options.addOption("type", true, "type of the GWAS data, 0 for simulated, 1 for PLINK tped format, default 0.");
        options.addOption("pathO", true, "path of the file recording the results, default 'result.txt''.");
        options.addOption("pathS2G", true, "the file containing mapping from SNPs to Genes, default NULL.");
        options.addOption("pathAG", true, "the file containing gene association, default NULL.");
        options.addOption("o", true, "the maximum order of SNP combinations in SHEIB, default -1.");
        options.addOption("maxGen", true, "the number of SNP combinations generated in SHEIB, default -1.");
        options.addOption("pb", true, "the probability of considering bioinformation in SHEIB, default 0.8.");
        options.addOption("cG", true, "threshold, default cG=0.05.");
        options.addOption("cGc", true, "threshold, default cGc=0.05.");
        options.addOption("rn", true, "threshold, default rn=-1.");
        options.addOption("seed", true, "seed of random, default 0");
        options.addOption("help", false, "help");
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            e.printStackTrace();
        }
        String arg=null;
        if(cmd.hasOption("help")){
            HelpFormatter formatter=new HelpFormatter();
            formatter.printHelp("sheib",options,true);
            return;
        }
        arg=cmd.getOptionValue("pathIn");
        if(arg!=null)
            this.pathIn=arg;
        arg=cmd.getOptionValue("pathO");
        if(arg!=null)
            this.pathO=arg;
        arg=cmd.getOptionValue("pathS2G");
        if(arg!=null)
            this.pathS2G=arg;
        arg=cmd.getOptionValue("pathAG");
        if(arg!=null)
            this.pathAG=arg;
        arg=cmd.getOptionValue("o");
        if(arg!=null)
            this.o=Integer.parseInt(arg);
        arg=cmd.getOptionValue("maxGen");
        if(arg!=null)
            this.maxGen=Integer.parseInt(arg);
        arg=cmd.getOptionValue("pb");
        if(arg!=null)
            this.pb=Double.parseDouble(arg);
        arg=cmd.getOptionValue("cG");
        if(arg!=null)
            this.cG=Double.parseDouble(arg);
        arg=cmd.getOptionValue("cGc");
        if(arg!=null)
            this.cGc=Double.parseDouble(arg);
        arg=cmd.getOptionValue("rn");
        if(arg!=null)
            this.rn=Integer.parseInt(arg);
        arg=cmd.getOptionValue("type");
        if(arg!=null)
            this.type=Integer.parseInt(arg);
        arg=cmd.getOptionValue("seed");
        if(arg!=null)
            this.seed=Integer.parseInt(arg);
        if(this.maxGen<0){
            System.out.println("maxGen must > 0");
            HelpFormatter formatter=new HelpFormatter();
            formatter.printHelp("sheib",options,true);
            return;
        }
        else {
            ready = true;
        }
    }
    public final boolean loadData(){
        if(this.type==0) {
            mem = new MemSimulated(this.pathIn);
        }
        else if(this.type==1){
            mem=new MemPLINK(this.pathIn+".tped",this.pathIn+".tfam");
        }
        if(mem.isReady()){
            if(this.pathS2G==null){
                return true;
            }
            else{
                snp2Gene = new SNP2Gene(this.pathS2G, mem);
                if (snp2Gene.isReady()) {
                    if(this.pathAG==null){
                        return true;
                    }
                    else {
                        gene2Gene = new Gene2Gene(this.pathAG, snp2Gene.getGene2Index());
                        if (gene2Gene.isReady()) {
                            return true;
                        }
                        else{
                            return false;
                        }
                    }
                }
                else{
                    return false;
                }
            }
        }
        else{
            return false;
        }
    }
    public final void run() {
        this.results=new TreeSet<>();
        if(this.o<=0){
            this.o=mem.recommendOrder();
        }
        this.random=new Random(this.seed);
        int g=0;
        double meanK2=0;
        while(g<this.maxGen){
            int[] x=generateSNPCom();
            double k2=mem.computeK2(x);
            meanK2=(meanK2*g+k2)/(g+1);
            if(k2<=meanK2){
                detectOnCom(x,k2);
            }
            g++;
            //System.out.printf("gen %d: result size=%d\n",g,results.size());
        }
    }
    private void detectOnCom(int[] x,double k2){
        if(x.length==1){
            return;
        }
        int[] xReduced=new int[x.length-1];
        double k2R=0;
        boolean successReduced=false;
        for(int i=0;i<x.length;i++){
            for(int j=0;j<i;j++){
                xReduced[j]=x[j];
            }
            for(int j=i;j<xReduced.length;j++){
                xReduced[j]=x[j+1];
            }
            k2R= mem.computeK2(xReduced);
            if(k2R<k2){
                successReduced=true;
                break;
            }
        }
        if(successReduced){
            detectOnCom(xReduced,k2R);
        }
        else{
            double g=mem.computeG(x);
            if(g<=this.cG) {
                results.add(new SNPComWithG(x, g));
            }
        }
        /*
        int l=o;
        while(l>1){
            int[] xx=new int[l-1];
            boolean reduced=false;
            for(int i=0;i<l;i++){
                for(int j=0;j<i;j++){
                    xx[j]=x[j];
                }
                for(int j=i;j<l-1;j++){
                    xx[j]=x[j+1];
                }
                double k2xx= mem.computeK2(xx);
                if(k2<k2xx){
                    x=xx;
                    k2=k2xx;
                    l=l-1;
                    reduced=true;
                    break;
                }
            }
            if(!reduced){
                double g=mem.computeG(x);
                if(g<=this.cG) {
                    results.add(new SNPCom(x, k2, g));
                }
            }
        }
         */
    }
    private final int[] generateSNPCom(){
        int[] r=new int[this.o];
        int ti=random.nextInt(mem.n);
        HashSet<Integer> snps=new HashSet<>();
        for(int i=0;i<this.o;i++){
            if(this.snp2Gene!=null&&this.snp2Gene.getGenesOfSNP(ti)!=null){
                for(int gene:this.snp2Gene.getGenesOfSNP(ti)){
                    snps.addAll(snp2Gene.getSNPsOfGene(gene));
                    if(gene2Gene!=null&&gene2Gene.getGenesOfGene(gene)!=null){
                        for(int gene2:gene2Gene.getGenesOfGene(gene)){
                            snps.addAll(snp2Gene.getSNPsOfGene(gene2));
                        }
                    }
                }
                for(int j=0;j<i;j++){
                    snps.remove(r[j]);
                }
            }
            if(snps.size()>0&&random.nextDouble()<this.pb){
                int randIndex=random.nextInt(snps.size());
                Iterator<Integer> it=snps.iterator();
                for(int j=0;j<randIndex;j++){
                    it.next();
                }
                ti=it.next();
                snps.remove(ti);
            }
            else{
                ti=random.nextInt(mem.n-i);
                for(int j=0;j<i;j++){
                    if(ti>=r[j]){
                        ti=ti+1;
                    }
                    else{
                        break;
                    }
                }
            }
            insertAsc(r,i,ti);
        }
        return r;
    }

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
    public void generateResults() throws IOException {
        BufferedOutputStream bos=new BufferedOutputStream(new FileOutputStream(this.pathO));
        Iterator<SNPComWithG> it=results.descendingIterator();
        while(it.hasNext()){
            bos.write(it.next().getOString(mem).getBytes());
        }
        bos.flush();
        bos.close();
    }
    public boolean isReady() {
        return ready;
    }
}
