package sly.gwas.algorithm;
import org.apache.commons.cli.*;
import sly.gwas.Funcs;
import sly.gwas.Logger;
import sly.gwas.mem.Mem;
import sly.gwas.mem.MemPLINK;
import sly.gwas.mem.MemSimulated;
import sly.gwas.snp.SNPComWithG;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * sparrow search algorithm：
 * 1 群智能算法；
 * 2 每代群体中有n只麻雀，每代选取种群中位置最好的pn只麻雀作为发现者，剩余的n-pn只麻雀作为跟随者；
 *   每次迭代中，每只麻雀都可以成为producer，根据适应度的值，但是，在整个种群中，producer的比例应该是不变的；
 * 3 发现者与追随者的更新公式不同：
 *   producers、scroungers
 *   发现者（producers）着重在它周围进行随机探查:
 *     随机生成一个R2，与安全阈值进行比较，以概率决定本代中发现者的更新方式：
 *       安全，着重在发现者的周围进行探索。
 *       不安全，着重在发现者的周围进行探索。
 *       看不懂，也觉得作者的思想有一些问题，所以，决定用字面意思来理解，安全近点，不安全，远点。
 *   更新群体的fitness和排序。
 *   追随者（explorer）着重跟着发现者跑：
 *     如果位置大于n/2，当前个体的更新方式着重于远离当前群体的最差个体。
 *     否则，当前个体的更新方式着重于靠近当前群体的最优个体。
 *   ？？选出群体中20%的个体
 *     如果比全局最优差，向全局最优靠近一点点。
 *     否则，根据当前个体与全局最差之间的距离，离全局最差远点。
 * 创新点：
 * 1 基于麻雀智能优化算法设计；
 * 2 SNP访问权重向量；
 * 3 针对上位性检测的特殊性，根据群体中基因的比例，决定采用的检测策略；
 * 4 解决假阳性问题。
 *
 */
public class SSA {
    class Spa{
        int[] x=null;
        double ce;
        double gini;
        double k2;
        int rankCe=0;
        int rankGini=0;
        int rankK2=0;
        int rankSum=0;
        public Spa(int [] x){
            if(!Funcs.check(x)){
                Logger.error("hehe\n");
            }
            //Funcs.testSolution(x);
            this.x=x;
            int[][] ca=mem.getTable(this.x);
            this.ce=mem.computeAdjustCe(ca,maxL);
            this.gini=mem.computeAdjustGini(ca,maxL);
            this.k2=mem.computeAdjustK2(ca,maxL);
        }
        public void remove(){
            spas.remove(this);
            for(Spa s:spas){
                if(this.ce<s.ce){
                    s.rankCe--;
                    s.rankSum--;
                }
                else if(this.ce>s.ce){
                    this.rankCe--;
                    this.rankSum--;
                }
                if(this.gini<s.gini){
                    s.rankGini--;
                    s.rankSum--;
                }
                else if(this.gini>s.gini){
                    this.rankGini--;
                    this.rankSum--;
                }
                if(this.k2<s.k2){
                    s.rankK2--;
                    s.rankSum--;
                }
                else if(this.k2>s.k2){
                    this.rankK2--;
                    this.rankSum--;
                }
            }
        }
        public void add2Spas(){
            for(Spa s:spas){
                if(this.ce<s.ce){
                    s.rankCe++;
                    s.rankSum++;
                }
                else if(this.ce>s.ce){
                    this.rankCe++;
                    this.rankSum++;
                }
                if(this.gini<s.gini){
                    s.rankGini++;
                    s.rankSum++;
                }
                else if(this.gini>s.gini){
                    this.rankGini++;
                    this.rankSum++;
                }
                if(this.k2<s.k2){
                    s.rankK2++;
                    s.rankSum++;
                }
                else if(this.k2>s.k2){
                    this.rankK2++;
                    this.rankSum++;
                }
            }
            for(int i=0;i<spas.size();i++){
                if(this.rankSum<spas.get(i).rankSum){
                    spas.add(i,this);
                    return;
                }
            }
            spas.add(this);
        }
        @Override
        public String toString() {
            return "Spa{" +
                    "x=" + Arrays.toString(x) +
                    ", ce=" + ce +
                    ", gini=" + gini +
                    ", k2=" + k2 +
                    ", rankCe=" + rankCe +
                    ", rankGini=" + rankGini +
                    ", rankK2=" + rankK2 +
                    ", rankSum=" + rankSum +
                    '}';
        }

        /**
         * 在本对象指代的SNP组合上检测上位性；
         */
        public void detectEpis(){
            int [] r=mem.detectEpisBasedOnMaxReducedK2(this.x);
            for(int i=0;i<d;i++){
                pss[this.x[i]]=pss[this.x[i]]*0.9;
            }
            if(r!=null) {
                SNPComWithG sg=new SNPComWithG(r,0);
                if(!results.contains(sg)) {
                    nTest++;
                    sg.setG(mem.computeG(r));
                    if (sg.getG() < cG) {
                        results.add(sg);
                    }
                }
            }
        }
    }
    /**
     * 算法参数；
     */
    int n;//n defines the number of sparrows
    int d=0;//the dimension of choice variables
    int maxG;//the maximum iterations
    double pd=0.4;//producer ratio (PD)
    double sd=0.2;//the ratio of sparrows who perceive the danger
    double st=0.8;//ST (safe threshold) is a value of 0.5 to 1.0
    int maxL=0;
    int seed=0;
    String pathIn="data.txt";
    String pathO="result.txt";
    int type=0;
    double cG=0.05;
    double thresholdSpasChaos=0.6;//越高，越倾向于变异，越低，越倾向于收敛。
    private int debugIndex=0;
    private int stepZero=2;
    /**
     * 算法运行的临时变量；
     */
    public boolean ready=false;
    ArrayList<Spa> spas=null;
    Mem mem=null;
    Random random=null;
    double[] pss=null;// the probability of each SNP being selected
    TreeSet<SNPComWithG> results=null;
    int nTest=0;
    @Override
    public String toString() {
        return "SSA{" +
                "n=" + n +
                ", d=" + d +
                ", maxG=" + maxG +
                ", pd=" + pd +
                ", sd=" + sd +
                ", st=" + st +
                ", maxL=" + maxL +
                ", seed=" + seed +
                ", pathIn='" + pathIn + '\'' +
                ", pathO='" + pathO + '\'' +
                ", type=" + type +
                ", cG=" + cG +
                '}';
    }

    void setParameters(String[] args){
        ready=false;
        Options options = new Options();
        options.addOption("help", false, "help");
        options.addOption("debugIndex", true, "debug GWAS index.");
        options.addOption("pathIn", true, "path of the GWAS data, default 'data.txt'.");
        options.addOption("pathO", true, "path of the file recording the results, default 'result.txt''.");
        options.addOption("type", true, "type of the GWAS data, 0 for simulated, 1 for PLINK tped format, default 0.");
        options.addOption("seed", true, "seed of random, default 0");
        options.addOption("d", true, "the dimension of SNP combinations in sparrows, default 0.");
        options.addOption("maxL", true, "the maximum length of contingency table, default 0.");
        options.addOption("cG", true, "threshold, default cG=0.05.");
        options.addOption("maxG", true, "the maximum iterations.");
        options.addOption("n", true, "the number of sparrows");
        options.addOption("pd", true, "the producer ratio (PD)");
        options.addOption("sd", true, "the ratio of sparrows who perceive the danger");
        options.addOption("r2", true, "R2 (alert value) is a number in the range of 0 to 1");
        options.addOption("st", true, "ST (safe threshold) is a value of 0.5 to 1.0");
        options.addOption("t", true, "thresholdSpasChaos, default 0.6.");

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
            formatter.printHelp("DE-SSA",options,true);
            return;
        }
        arg=cmd.getOptionValue("debugIndex");
        if(arg!=null)
            this.debugIndex=Integer.parseInt(arg);
        arg=cmd.getOptionValue("n");
        if(arg!=null)
            this.n=Integer.parseInt(arg);
        arg=cmd.getOptionValue("d");
        if(arg!=null)
            this.d=Integer.parseInt(arg);
        arg=cmd.getOptionValue("maxG");
        if(arg!=null)
            this.maxG=Integer.parseInt(arg);
        arg=cmd.getOptionValue("pd");
        if(arg!=null)
            this.pd=Double.parseDouble(arg);
        arg=cmd.getOptionValue("sd");
        if(arg!=null)
            this.sd=Double.parseDouble(arg);
        arg=cmd.getOptionValue("st");
        if(arg!=null)
            this.st=Double.parseDouble(arg);
        arg=cmd.getOptionValue("maxL");
        if(arg!=null)
            this.maxL=Integer.parseInt(arg);
        arg=cmd.getOptionValue("seed");
        if(arg!=null)
            this.seed=Integer.parseInt(arg);
        arg=cmd.getOptionValue("pathIn");
        if(arg!=null)
            this.pathIn=arg;
        arg=cmd.getOptionValue("pathO");
        if(arg!=null)
            this.pathO=arg;
        arg=cmd.getOptionValue("type");
        if(arg!=null)
            this.type=Integer.parseInt(arg);
        arg=cmd.getOptionValue("cG");
        if(arg!=null)
            this.cG=Double.parseDouble(arg);
        arg=cmd.getOptionValue("t");
        if(arg!=null)
            this.thresholdSpasChaos=Double.parseDouble(arg);
        if(this.maxG<0){
            Logger.error("maxGen must > 0");
            HelpFormatter formatter=new HelpFormatter();
            formatter.printHelp("DE-SSA",options,true);
            return;
        }
        else {
            ready = true;
        }

    }

    public SSA(String[] args){
        this.setParameters(args);
        this.random=new Random(this.seed);
    }

    /**
     * 初始化麻雀群体的位置；
     */
    void init(){
        this.spas=new ArrayList<>(n);
        for(int i=0;i<this.n;i++){
            int[] x=genRanSpa();
            Spa spa=new Spa(x);
            spa.add2Spas();
        }
        this.pss=new double[this.mem.n];
        this.results=new TreeSet<>();
        Arrays.fill(this.pss,1.0);
        nTest=0;
    }

    /**
     * 随机生成一只麻雀；
     */
    private int[] genRanSpa(){
        int[] r=Funcs.sample(random,this.mem.n,this.d);
        return r;
    }

    /**
     * 从外存载入数据；
     * @return 是否成功；
     */
    public final boolean loadData() {
        if (this.type == 0) {
            mem = new MemSimulated(this.pathIn);
        } else if (this.type == 1) {
            mem = new MemPLINK(this.pathIn + ".tped", this.pathIn + ".tfam");
        }
        if (mem.isReady()) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * SSA算法的主要流程控制代码；
     */
    public void run(){
        //fix d,maxL
        if(this.d<=0){
            this.d=mem.recommendOrder();
        }
        if(this.maxL<=0){
            this.maxL=mem.recommendMaxL();
        }
        this.init();
        //detectEpis();
        for(int i=0;i<maxG;i++){
            this.updateSSAPositions();
            HashSet<Integer> snps=new HashSet<>();
            for(Spa spa:spas){
                for(int s:spa.x){
                    snps.add(s);
                }
            }
            if(snps.size()<n*thresholdSpasChaos*d){
                detectEpis(true);
            }
            else{
                detectEpis(false);
            }
            if(debugIndex>0&&i%debugIndex==0){
                this.generateDebugResults(i);
            }
        }
    }
    private final void detectEpis(boolean removeHead){
        if(removeHead) {
            for (int i = 0; i < n * this.sd; i++) {
                Spa spa = this.spas.get(0);
                spa.detectEpis();
                spa.remove();
            }
        }
        else{
            for (int i = 0; i < n * this.sd; i++) {
                Spa spa = this.spas.get(i);
                spa.detectEpis();
            }
            for(int i=0;i<n*this.sd;i++){
                Spa spa = this.spas.get(this.spas.size()-1);
                spa.remove();
            }
        }
    }
    /**
     * 更新麻雀群体中每个个体的位置；
     */
    private void updateSSAPositions(){
        Logger.debug("updateSSAPositions\n");
        int pd=(int)(this.pd*this.n);
        /**
         * 发现者更新方法：
         * 随机生成一个随机数，以概率ST选择两种移动方式之一：
         * 1 模拟麻雀在附近移动：随机选择个体中的一个维度的值，进行随机替换；
         * 2 模拟麻雀在附近移动的幅度比较大：随机选择个体一半的维度，进行随机替换；
         */
        if(random.nextDouble()<this.st){
            for(int i=0;i<pd;i++){
                int[] x=spas.get(i).x;
                int[] nx=new int[this.d];
                int indexReplace=random.nextInt(this.d);
                for(int j=0;j<indexReplace;j++){
                    nx[j]=x[j];
                }
                for(int j=indexReplace+1;j<this.d;j++){
                    nx[j-1]=x[j];
                }
                Funcs.insertAsc(this.random,nx,this.mem.n,this.d-1);
                Spa s=new Spa(nx);
                Spa so=spas.get(i);
                s.add2Spas();
                if(s.rankSum<so.rankSum){
                    so.remove();
                }
                else{
                    s.remove();
                }
            }
        }
        else{
            for(int i=0;i<pd;i++){
                int[] x=spas.get(i).x;
                int[] nx=new int[this.d];
                int nReplace=this.d/2;
                int jx=0;
                int jnx=0;
                for(int indexReplace:Funcs.sample(random,this.d,nReplace)){
                    while(jx<indexReplace){
                        nx[jnx++]=x[jx++];
                    }
                    jx++;
                }
                while(jx<this.d){
                    nx[jnx++]=x[jx++];
                }
                for(int j=0;j<nReplace;j++) {
                    Funcs.insertAsc(this.random, nx, this.mem.n,  jnx++);
                }
                Spa s=new Spa(nx);
                Spa so=spas.get(i);
                s.add2Spas();
                if(s.rankSum<so.rankSum){
                    so.remove();
                }
                else{
                    s.remove();
                }
            }
        }
        /**
         * 追随者更新方式，根据追随者的位置序号是否在前n/2，决定不同的移动方式：
         * 1 前n/2的麻雀，更接近中心，所以它们向最优的个体移动；
         * 2 后n/2的麻雀，远离群体中的最差值；
         */
        for(int i=pd;i<n;i++){
            if(i<n/2){
                int selectProducer=Math.min(random.nextInt(pd),random.nextInt(pd));
                int [] nx=Funcs.add(random,spas.get(i).x,spas.get(selectProducer).x);
                Spa s=new Spa(nx);
                Spa so=spas.get(i);
                s.add2Spas();
                if(s.rankSum<so.rankSum){
                    so.remove();
                }
                else{
                    s.remove();
                }
            }
            else{
                int [] nx=Funcs.sample(random,this.mem.n,this.d,this.pss);
                Spa s=new Spa(nx);
                Spa so=spas.get(i);
                s.add2Spas();
                if(s.rankSum<so.rankSum){
                    so.remove();
                }
                else{
                    s.remove();
                }
            }
        }
        /**
         * 从群体中随机选择一定比例的个体做如下操作：
         * 1 如果选择的是最优个体，最优个体进行重新生成；
         * 2 如果选择的不是最优个体，将这个个体向最优个体进行靠近。
         */
        for(int i:Funcs.sample(random,n,(int)(n*this.sd))){
            if(i==0){
                int [] nx=Funcs.sample(random,this.mem.n,this.d,this.pss);
                //spas.get(i).remove();
                Spa s=new Spa(nx);
                s.add2Spas();
            }
            else{
                int t=Math.min(pd,i+1);
                int selectProducer=Math.min(random.nextInt(t),random.nextInt(t));
                int [] nx=Funcs.add(random,spas.get(i).x,spas.get(selectProducer).x);
                //spas.get(i).remove();
                Spa s=new Spa(nx);
                s.add2Spas();
            }
        }
    }

    /**
     * 将检测到的上位性存储到输出文件中；
     * 利用四分位差降低假阳性；
     * @throws IOException 文件IO异常。
     */
    public void generateResults2() throws IOException {
        if(results.size()>=4){
            int q3Index=results.size()/4*3;
            int q1Index=results.size()/4;
            double q1=0;
            double q3=0;
            Iterator<SNPComWithG> it=results.descendingIterator();
            for(int i=0;it.hasNext();i++){
                if(i==q1Index){
                    q1=it.next().getG();
                }
                else if(i==q3Index){
                    q3=it.next().getG();
                }
                else{
                    it.next();
                }
            }
            double iqr=q3-q1;
            double bound=q1-1.5*iqr;
            it=results.descendingIterator();
            BufferedOutputStream bos=new BufferedOutputStream(new FileOutputStream(this.pathO));
            while(it.hasNext()){
                SNPComWithG s=it.next();
                if(s.getG()<bound){
                    bos.write(s.getOString(mem).getBytes());
                }
            }
            bos.flush();
            bos.close();
        }
    }

    /**
     * 将检测到的上位性存储到输出文件中；
     * 利用数字上的跨度异常来降低假阳性；
     * @throws IOException 文件IO异常。
     */
    public void generateResults() throws IOException {
       if(results.size()<=1){
           Iterator<SNPComWithG> it=results.descendingIterator();
           BufferedOutputStream bos=new BufferedOutputStream(new FileOutputStream(this.pathO));
           while(it.hasNext()){
               bos.write(it.next().getOString(mem).getBytes());
           }
           bos.flush();
           bos.close();
       }
       else{
           double biggestStep=0;
           int biggestStepIndex=0;
           Iterator<SNPComWithG> it=results.descendingIterator();
           double g1=it.next().getG();
           int i1=0;
           while(it.hasNext()){

               double g2=it.next().getG();
               double step=Math.pow(10,this.stepZero);
               if(g1!=0){
                   step=g2/g1;
               }
               if(step>=biggestStep){
                   biggestStep=step;
                   biggestStepIndex=i1;
               }
               g1=g2;
               i1++;
               /*
               double g2=it.next().getG();
               double step=1;
               if(g1!=0){
                   step=g2/g1;
               }
               if(step>=biggestStep){
                   biggestStep=step;
                   biggestStepIndex=i1;
               }
               g1=g2;
               i1++;
               */
           }
           it=results.descendingIterator();
           BufferedOutputStream bos=new BufferedOutputStream(new FileOutputStream(this.pathO));
           for(int i=0;i<=biggestStepIndex;i++){
               bos.write(it.next().getOString(mem).getBytes());
           }
           bos.flush();
           bos.close();
       }
    }

    /**
     * 针对算法在真实数据上的检测，提供调试结果：
     *   起因在于，在真实数据集上算法运行的时间太长了，无法保证服务器可以一次完整运行完毕；
     *   因此算法的迭代过程中，会每迭代debugIndex次，进行一次结果输出，并将算法的当前结果输出到debug目录下；
     * @param i 当前迭代次数。
     */
    private void generateDebugResults(int i){
        File file=new File("debug");
        if(!file.exists()){
            file.mkdir();
        }
        file=new File(String.format("debug/%d",i));
        if(!file.exists()){
            file.mkdir();
        }
        String pathO=this.pathO;
        this.pathO=String.format("debug/%d/%s",i,new File(this.pathIn).getName());
        try {
            generateResults();
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.pathO=pathO;
    }
}
