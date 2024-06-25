package sly.gwas.algorithm;
import org.apache.commons.cli.*;
import sly.gwas.Logger;
import sly.gwas.mem.Mem;
import sly.gwas.mem.MemPLINK;
import sly.gwas.mem.MemSimulated;
import sly.gwas.snp.SNPCom;
import sly.gwas.snp.SNPComWithG;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * 灰狼算法在GWAS上的实现；
 */
public class GWO {
    String pathIn="data.txt";
    String pathO="result.txt";
    int o=-1;
    int maxL=-1;
    int numWolves=40;
    int[][] wolves=null;
    double[] leaderWeights=null;
    int maxGen=-1;
    double cG=0.05;
    int type=0;
    int seed=0;
    //头狼一共有3只，记录狼群中不同评测指标下最优秀的狼，依次为:ace,agini,ak2；
    int leaderACe=0;
    int leaderAGini=0;
    int leaderAK2=0;
    double[] aCes=null;
    double[] aGinis=null;
    double[] aK2s=null;
    Mem mem=null;
    Random random=null;
    boolean ready=false;
    TreeSet<SNPComWithG> results=null;
    int gen=0;
    TreeSet<SNPCom> history=null;
    TreeMap<SNPCom,Integer> historyLeaders=null;
    double maxRate=1;
    double varRate=0;
    double leaderIncreaseWeight=0.2;
    boolean debugGwas=false;

    public int getLl() {
        return ll;
    }

    public void setLl(int ll) {
        this.ll = ll;
    }

    //日志级别
    int ll=0;

    /**
     * 灰狼算法构造函数；
     * @param args 参数列表；
     */
    public GWO(String [] args){
        ready=false;
        Options options = new Options();
        options.addOption("help", false, "help");
        options.addOption("pathIn", true, "path of the GWAS data, default 'data.txt'.");
        options.addOption("pathO", true, "path of the file recording the results, default 'result.txt''.");
        options.addOption("type", true, "type of the GWAS data, 0 for simulated, 1 for PLINK tped format, default 0.");
        options.addOption("seed", true, "seed of random, default 0");
        options.addOption("o", true, "the maximum order of SNP combinations, default -1.");
        options.addOption("maxL", true, "the maximum length of contingency table, default -1.");
        options.addOption("cG", true, "threshold, default cG=0.05.");
        options.addOption("maxGen", true, "the number of SNP combinations generated, default -1.");
        options.addOption("numWolves", true, "number of the wolves in the algorithm");
        options.addOption("leaderIncreaseWeight", true, "leaderIncreaseWeight, default 0.2");
        options.addOption("maxRate", true, "maxRate, default 1");
        options.addOption("ll", true, "logger level, 0 for info, 1 for debug");
        options.addOption("debugGwas", false, "help");

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
            formatter.printHelp("de-gwo",options,true);
            return;
        }
        if(cmd.hasOption("debugGwas")){
            this.debugGwas=true;
        }
        arg=cmd.getOptionValue("pathIn");
        if(arg!=null)
            this.pathIn=arg;
        arg=cmd.getOptionValue("pathO");
        if(arg!=null)
            this.pathO=arg;
        arg=cmd.getOptionValue("type");
        if(arg!=null)
            this.type=Integer.parseInt(arg);
        arg=cmd.getOptionValue("seed");
        if(arg!=null)
            this.seed=Integer.parseInt(arg);
        arg=cmd.getOptionValue("o");
        if(arg!=null)
            this.o=Integer.parseInt(arg);
        arg=cmd.getOptionValue("maxL");
        if(arg!=null)
            this.maxL=Integer.parseInt(arg);
        arg=cmd.getOptionValue("cG");
        if(arg!=null)
            this.cG=Double.parseDouble(arg);
        arg=cmd.getOptionValue("maxGen");
        if(arg!=null)
            this.maxGen=Integer.parseInt(arg);
        arg=cmd.getOptionValue("numWolves");
        if(arg!=null)
            this.numWolves=Integer.parseInt(arg);
        arg=cmd.getOptionValue("leaderIncreaseWeight");
        if(arg!=null)
            this.leaderIncreaseWeight=Double.parseDouble(arg);
        arg=cmd.getOptionValue("maxRate");
        if(arg!=null)
            this.maxRate=Double.parseDouble(arg);
        arg=cmd.getOptionValue("ll");
        if(arg!=null)
            this.ll=Integer.parseInt(arg);
        if(this.maxGen<0){
            Logger.error("maxGen must > 0");
            HelpFormatter formatter=new HelpFormatter();
            formatter.printHelp("DE-GWO",options,true);
            return;
        }
        else {
            ready = true;
        }

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
     * 算法执行主函数；
     */
    public final void run() {
        this.results=new TreeSet<>();
        if(this.o<=0){
            this.o=mem.recommendOrder();
        }
        if(this.maxL<=0){
            this.maxL=mem.recommendMaxL();
        }

        StringBuffer sb=new StringBuffer();
        sb.append("DE-GWO(");
        sb.append("pathIn: ");
        sb.append(pathIn);
        sb.append(",");
        sb.append("pathO: ");
        sb.append(pathO);
        sb.append(",");
        sb.append("type: ");
        sb.append(type);
        sb.append(",");
        sb.append("seed: ");
        sb.append(seed);
        sb.append(",");
        sb.append("o: ");
        sb.append(o);
        sb.append(",");
        sb.append("maxL: ");
        sb.append(maxL);
        sb.append(",");
        sb.append("cG: ");
        sb.append(cG);
        sb.append(",");
        sb.append("maxGen: ");
        sb.append(maxGen);
        sb.append(",");
        sb.append("numWolves: ");
        sb.append(numWolves);
        sb.append(",");
        sb.append("leaderIncreaseWeight: ");
        sb.append(leaderIncreaseWeight);
        sb.append(",");
        sb.append("maxRate: ");
        sb.append(maxRate);
        sb.append(")\n");
        Logger.info(sb.toString());
        this.random=new Random(this.seed);
        this.gen=0;
        //初始化狼群，并计算适应度；
        initWolves();
        updateVarRate();
        history=new TreeSet<>();
        historyLeaders=new TreeMap<>();
        while(true){
            //选择头狼；
            selectLeaders();
            print();
            //在头狼上检测上位性；
            detectOnCom(wolves[leaderACe],mem.computeK2(wolves[leaderACe]));
            detectOnCom(wolves[leaderAGini],mem.computeK2(wolves[leaderAGini]));
            detectOnCom(wolves[leaderAK2],mem.computeK2(wolves[leaderAK2]));
            this.gen++;
            if(this.gen==this.maxGen){
                break;
            }
            if(this.debugGwas){
                if(this.gen%20000==0){
                    this.generateDebugResults();
                }
            }
            //向头狼移动；
            moveWolves2Leaders();
            //System.out.printf("gen %d: result size=%d\n",g,results.size());
            //更新变异率；
            updateVarRate();
        }
    }

    /**
     * 更新变异率；
     * 变异率和当前狼群中的复杂度挂钩；
     */
    private void updateVarRate(){
        HashSet<Integer> species=new HashSet<>();
        for(int[] x:wolves){
            for(int i:x){
                species.add(i);
            }
        }
        varRate=(double)species.size();
        varRate=(1-(varRate-o)/((numWolves-1)*o))*maxRate;
    }
    /**
     * 将狼群向头狼进行随机移动；
     */
    private void moveWolves2Leaders(){
        for(int i=0;i<this.numWolves;i++){
            if(!(i==this.leaderACe||i==this.leaderAGini||i==this.leaderAK2)){
                moveWolf2Leaders(i);
            }
        }
    }

    /**
     * 将一头狼向头狼进行移动；
     * @param wolf 狼下标；
     */
    private void moveWolf2Leaders(int wolf){
        //System.out.println(String.format("update wolves %d",wolf));
        int[] x1=move(wolves[wolf],wolves[leaderACe]);
        int[] x2=move(wolves[wolf],wolves[leaderAGini]);
        int[] x3=move(wolves[wolf],wolves[leaderAK2]);
        wolves[wolf]=meanWolves(x1,x2,x3);
        int [][] ca=mem.getTable(wolves[wolf]);
        aCes[wolf]=mem.computeAdjustCe(ca,maxL);
        aGinis[wolf]=mem.computeAdjustGini(ca,maxL);
        aK2s[wolf]=mem.computeAdjustK2(ca,maxL);
        SNPCom s=new SNPCom(wolves[wolf]);
        /*
        if(s.contains(new int[]{0,49,97})){
            Logger.debug("##################\n");
            Logger.debug(Arrays.toString(wolves[wolf]));
        }

         */
        if(historyLeaders.containsKey(s)){
            leaderWeights[wolf]=1+leaderIncreaseWeight*historyLeaders.get(s);
        }
        else{
            leaderWeights[wolf]=1;
        }
    }

    /**
     * 取三匹狼位置的均值；
     * @param a 狼；
     * @param b 狼；
     * @param c 狼；
     * @return 新的位置；
     */
    private int[] meanWolves(int[] a,int[] b,int[] c){
        int[] r=new int[this.o];
        //记录a中被抽取的元素的地址；
        boolean[] selectedInA=new boolean[this.o];
        Arrays.fill(selectedInA,false);
        int countInA=0;
        //记录b中被抽取的元素的地址；
        boolean[] selectedInB=new boolean[this.o];
        Arrays.fill(selectedInB,false);
        int countInB=0;
        //记录c中被抽取的元素的地址；
        boolean[] selectedInC=new boolean[this.o];
        Arrays.fill(selectedInC,false);
        int countInC=0;
        //从a中选择1/3的元素
        for(int i=0;i<this.o/3;i++){
            int t=random.nextInt(this.o-countInA);
            while(selectedInA[t]){
                t++;
            }
            int selectedValue=a[t];
            countInA++;
            selectedInA[t]=true;
            for(int j=0;j<this.o;j++){
                if(b[j]==selectedValue){
                    selectedInB[j]=true;
                    countInB++;
                    break;
                }
            }
            for(int j=0;j<this.o;j++){
                if(c[j]==selectedValue){
                    selectedInC[j]=true;
                    countInC++;
                    break;
                }
            }
            insertAsc(r,i,selectedValue);
        }
        //从b中选择1/3的元素
        for(int i=this.o/3;i<this.o*2/3;i++){
            int t=random.nextInt(this.o-countInB);
            while(selectedInB[t]){
                t++;
            }
            int selectedValue=b[t];
            countInB++;
            selectedInB[t]=true;
            for(int j=0;j<this.o;j++){
                if(c[j]==selectedValue){
                    selectedInC[j]=true;
                    countInC++;
                    break;
                }
            }
            insertAsc(r,i,selectedValue);
        }
        //从c中选择1/3的元素
        for(int i=this.o*2/3;i<this.o;i++){
            int t=random.nextInt(this.o-countInC);
            while(selectedInC[t]){
                t++;
            }
            int selectedValue=c[t];
            countInC++;
            selectedInC[t]=true;
            insertAsc(r,i,selectedValue);
        }
        return r;
    }

    /**
     * 将狼a向狼b移动，生成新的位置；
     * @param a 狼；
     * @param b 狼；
     * @return 新的位置；
     */
    private int[] move(int[] a,int[] b){
        int[] r=new int[this.o];
        //记录a中被抽取的元素的地址；
        boolean[] selectedInA=new boolean[this.o];
        Arrays.fill(selectedInA,false);
        int countInA=0;
        //记录b中被抽取的元素的地址；
        boolean[] selectedInB=new boolean[this.o];
        Arrays.fill(selectedInB,false);
        int countInB=0;
        //从a中选择一半的元素
        for(int i=0;i<this.o/2;i++){
            if(random.nextDouble()<=varRate){
                int selectedValue=random.nextInt(mem.n-i);
                for(int j=0;j<i;j++){
                    if(selectedValue>=r[j]){
                        selectedValue=selectedValue+1;
                    }
                    else{
                        break;
                    }
                }for (int j = 0; j < this.o; j++) {
                    if (a[j] == selectedValue) {
                        selectedInA[j] = true;
                        countInA++;
                        break;
                    }
                }
                for (int j = 0; j < this.o; j++) {
                    if (b[j] == selectedValue) {
                        selectedInB[j] = true;
                        countInB++;
                        break;
                    }
                }
                insertAsc(r,i,selectedValue);
            }
            else {
                int t = random.nextInt(this.o - countInA);
                while (selectedInA[t]) {
                    t++;
                }
                int selectedValue = a[t];
                countInA++;
                selectedInA[t] = true;
                for (int j = 0; j < this.o; j++) {
                    if (b[j] == selectedValue) {
                        selectedInB[j] = true;
                        countInB++;
                        break;
                    }
                }
                insertAsc(r, i, selectedValue);
            }
        }
        //从b中选择一半的元素
        for(int i=this.o/2;i<this.o;i++){
            if(random.nextDouble()<=varRate) {
                int selectedValue=random.nextInt(mem.n-i);
                for(int j=0;j<i;j++){
                    if(selectedValue>=r[j]){
                        selectedValue=selectedValue+1;
                    }
                    else{
                        break;
                    }
                }
                for (int j = 0; j < this.o; j++) {
                    if (b[j] == selectedValue) {
                        selectedInB[j] = true;
                        countInB++;
                        break;
                    }
                }
                insertAsc(r,i,selectedValue);
            }
            else{
                int t = random.nextInt(this.o - countInB);
                while (selectedInB[t]) {
                    t++;
                }
                int selectedValue = b[t];
                countInB++;
                selectedInB[t] = true;
                insertAsc(r, i, selectedValue);
            }
        }
        return r;
    }

    /**
     * 在SNP组合上检测上位性；
     * @param x SNP组合；
     * @param k2 x的k2值；
     */
    private void detectOnCom(int[] x,double k2){
        SNPCom s=new SNPCom(x);
        if(x.length==1||history.contains(s)){
            //System.out.println("his");
            return;
        }
        Logger.debug(String.format("detect on %s\n",Arrays.toString(x)));
        history.add(s);
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
                SNPComWithG snp=new SNPComWithG(x, g);
                results.add(snp);
                Logger.debug(snp.getOString(mem));
            }
        }
    }
    /**
     * 随机初始化狼群；
     * 并计算适应度；
     * 选择头狼；
     */
    private void initWolves(){
        this.wolves=new int[numWolves][];
        this.aCes=new double[this.numWolves];
        this.aGinis=new double[this.numWolves];
        this.aK2s=new double[this.numWolves];
        this.leaderWeights=new double[this.numWolves];
        Arrays.fill(this.leaderWeights,1);
        for(int i=0;i<numWolves;i++){
            this.wolves[i]=genRanWolf();
            int[][] ca=mem.getTable(this.wolves[i]);
            this.aCes[i]=mem.computeAdjustCe(ca,maxL);
            this.aGinis[i]=mem.computeAdjustGini(ca,maxL);
            this.aK2s[i]=mem.computeAdjustK2(ca,maxL);
        }
    }
    /**
     * 选举头狼；
     * 根据ce，gini和k2的值，分别在群体中选择最优秀的三个不同的个体作为头狼；
     * 在选择头狼时，需要考虑头狼的leaderweights，确保乘以该值的情况下，选择的头狼最优。
     */
    private void selectLeaders(){
        leaderACe=0;
        for(int i=1;i<numWolves;i++){
            if(aCes[i]*leaderWeights[i]<aCes[leaderACe]*leaderWeights[leaderACe]){
                leaderACe=i;
            }
        }
        leaderAGini=0;
        while(leaderAGini==leaderACe){
            leaderAGini++;
        }
        for(int i=1;i<numWolves;i++){
            if(i!=leaderACe&&aGinis[i]*leaderWeights[i]<aGinis[leaderAGini]*leaderWeights[leaderAGini]){
                leaderAGini=i;
            }
        }
        leaderAK2=0;
        while(leaderAK2==leaderACe||leaderAK2==leaderAGini){
            leaderAK2++;
        }
        for(int i=1;i<numWolves;i++){
            if(i!=leaderACe&&i!=leaderAGini&&aK2s[i]*leaderWeights[i]<aK2s[leaderAK2]*leaderWeights[leaderAK2]){
                leaderAK2=i;
            }
        }
        /*
        SNPCom s=new SNPCom(this.wolves[leaderACe]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderACe]=1+leaderIncreaseWeight*historyLeaders.get(s);
        s=new SNPCom(this.wolves[leaderAGini]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderACe]=1+leaderIncreaseWeight*historyLeaders.get(s);
        s=new SNPCom(this.wolves[leaderAK2]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderACe]=1+leaderIncreaseWeight*historyLeaders.get(s);
        */
        SNPCom s=new SNPCom(this.wolves[leaderACe]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderACe]=1+leaderIncreaseWeight*historyLeaders.get(s);
        s=new SNPCom(this.wolves[leaderAGini]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderAGini]=1+leaderIncreaseWeight*historyLeaders.get(s);
        s=new SNPCom(this.wolves[leaderAK2]);
        historyLeaders.put(s,historyLeaders.getOrDefault(s,0)+1);
        this.leaderWeights[leaderAK2]=1+leaderIncreaseWeight*historyLeaders.get(s);
    }
    /**
     * 随机生成一匹狼；
     */
    private int[] genRanWolf(){
        int[] r=new int[this.o];
        int ti=-1;
        for(int i=0;i<this.o;i++){
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
     * 打印当前群体信息；
     */
    public void print(){
        if(this.ll==0){
            if(this.gen%2000==0) {
                StringBuffer sb = new StringBuffer();
                sb.append("gen ");
                sb.append(this.gen);
                sb.append(" result size ");
                sb.append(this.results.size());
                sb.append("\n");
                Logger.info(sb.toString());
            }
        }
        else {
            StringBuffer sb=new StringBuffer();
            sb.append("gen ");
            sb.append(this.gen);
            sb.append(":\n");
            for (int i : new int[]{leaderACe, leaderAGini, leaderAK2}) {
                sb.append(i);
                sb.append(Arrays.toString(wolves[i]));
                sb.append(" ");
                sb.append(aCes[i]);
                sb.append(",");
                sb.append(aGinis[i]);
                sb.append(",");
                sb.append(aK2s[i]);
                sb.append("\n");
            }
            sb.append("varRate: ");
            sb.append(varRate);
            sb.append(" result size: ");
            sb.append(this.results.size());
            sb.append("\n");
            Logger.debug(sb.toString());
        }
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
    private void generateDebugResults(){
        File dir=new File("debug-gwas");
        if (!dir.exists()){
            dir.mkdir();
        }
        dir=new File(String.format("debug-gwas/%d",this.gen));
        if (!dir.exists()){
            dir.mkdir();
        }
        try {
            BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(String.format("debug-gwas/%d/%s", this.gen, new File(this.pathO).getName())));
            Iterator<SNPComWithG> it = results.descendingIterator();
            while (it.hasNext()) {
                bos.write(it.next().getOString(mem).getBytes());
            }
            bos.flush();
            bos.close();
        }
        catch (IOException e){
            Logger.error("generateDebugResults error");
        }
    }
    public boolean isReady() {
        return ready;
    }
}
