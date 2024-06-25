package sly.gwas;
import sly.gwas.algorithm.GWO;
import sly.gwas.algorithm.SSA;

import java.io.IOException;
import java.rmi.server.ExportException;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/bd_gwas -pathO bd_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/cad_gwas -pathO cad_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/cd_gwas -pathO cd_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/ht_gwas -pathO ht_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/ra_gwas -pathO ra_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/t1d_gwas -pathO t1d_gwas.txt &
 * java -jar /home/sly/gwas/baselines/EPISSA/method.jar -type 1 -maxG 6000000 -n 600 -debugIndex 20000 -pathIn /home/sly/gwas/data/wtccc/t2d_gwas -pathO t2d_gwas.txt &
 *
 */
public class Main {
    private static void mainGWO(String[] args){
        //String pathIn="E:\\data\\gwas\\sim_data_20220613\\DME3_1600_1000\\146_915.txt";
        //String pathIn="d:\\workspace\\java\\de-gwo\\example";
        //GWO gwo=new GWO(new String[]{"-pathIn",pathIn,"-maxGen","200000","-numWolves","20","-seed","0","-type","1","-ll","0"});
        GWO gwo=new GWO(args);
        switch (gwo.getLl()){
            case 0:
                Logger.setLevel(Logger.INFO);
                break;
            case 1:
                Logger.setLevel(Logger.DEBUG);
                break;
            default:
                Logger.setLevel(Logger.INFO);
                break;
        }
        if(gwo.isReady()) {
            if (gwo.loadData()) {
                gwo.run();
                try {
                    gwo.generateResults();
                } catch (IOException e) {
                    e.printStackTrace();
                    Logger.error("生成结果文件失败");
                }
            }
        }
    }
    private static void mainSSA(String [] args){
        //String pathIn="D:\\workspace\\gwas\\data\\sim-data-20230311\\DNME04_1600_1000\\173_339.txt";
        //String pathIn="example";
/*
        args=new String[]{
                "-pathIn",pathIn,
                "-maxG","8000",
                "-n","40",
                "-seed","0",
                "-type","1",
                "-pd","0.4",
                "-sd","0.2",
                "-st","0.8",
                "-debugIndex","40"
        };
*/
        SSA ssa=new SSA(args);
        if(ssa.ready){
            System.out.println(ssa);
            if(ssa.loadData()){
                ssa.run();
            }
            try{
                ssa.generateResults();
            }
            catch (Exception e){
                System.out.println(e);
            }
        }

    }
    public static void main(String[] args) {
        mainSSA(args);
    }
}
