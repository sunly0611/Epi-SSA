package sly.gwas.mem;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.io.FileUtils;
import sly.gwas.Logger;

/**
 * 存储维护模拟GWAS数据的类；
 */
final public class MemSimulated extends Mem {
    /**
     * 构造函数，根据文件路径，载入GWAS数据；
     * 模拟数据；
     * @param filename
     */
    public MemSimulated(String filename){
        super();
        Logger.info(String.format("MemSimulated loading %s\n", filename));
        File f=new File(filename);
        if(!f.exists()){
            Logger.error("file does not exist\n");
        }
        else{
            try {
                if (!fillMN(f)) {
                    Logger.error("file format error\n");
                } else {
                    Logger.info(String.format("n:%d case:%d control:%d\n", n, m1, m0));
                    if (!fillMem(f)) {
                        Logger.error("fillMem error\n");
                    } else {
                        ready = true;
                    }
                }
            }
            catch (Exception e){
                Logger.error("file load error\n");
            }
        }
    }
    /**
     * 根据GWAS数据文件的内容，获取，样本数目和SNP的数目；
     * @param file
     * @return
     */
    private final boolean fillMN(File file) throws IOException {
        List<String> lines = FileUtils.readLines(file, "UTF-8");
    	m0=0;
    	m1=0;
    	Iterator<String> it=lines.iterator();
    	n=it.next().split("\t").length-1;
    	while(it.hasNext()) {
    		String line=it.next();
    		String[] ss=line.split("\t");
    		if(ss[ss.length-1].equals("1")) {
    			m1++;
    		}
    		else if(ss[ss.length-1].equals("0")) {
    			m0++;
    		}
    		else {
    			System.out.printf("phenotype should be 0 or 1\n");
    			return false;
    		}
    	}
    	m=m0+m1;
        return true;
    }
    /**
     * 根据文件的内容，载入内存；
     * @param file GWAS模拟数据路径；
     * @return 是否成功载入数据；
     */
    private final boolean fillMem(File file) throws IOException {
    	l1=(int)Math.ceil(((double)m1)/SIZE);
        l0=(int)Math.ceil(((double)m0)/SIZE);
        mem1=new long[n][3][l1];
        mem0=new long[n][3][l0];
        names=new String[n];
        List<String> lines = FileUtils.readLines(file, "UTF-8");
        Iterator<String> it=lines.iterator();
        String[] ss=it.next().split("\t");
        System.arraycopy(ss,0,names,0,n);
        for(int i=0,i0=0,i1=0;i<m;i++){
        	ss=it.next().split("\t");
            if(ss[ss.length-1].equals("1")){
                int indexOfVector=i1/SIZE;
                int indexInVector=i1-SIZE*indexOfVector;
                i1++;
                for(int j=0;j<n;j++){
                    mem1[j][Integer.parseInt(ss[j])][indexOfVector]|=ONE>>>indexInVector;
                }
            }
            else if(ss[ss.length-1].equals("0")){
                int indexOfVector=i0/SIZE;
                int indexInVector=i0-SIZE*indexOfVector;
                i0++;
                for(int j=0;j<n;j++){
                	mem0[j][Integer.parseInt(ss[j])][indexOfVector]|=ONE>>>indexInVector;
                }
            }
            else{
                System.out.printf("phenotype should be 0 or 1\n");
                return false;
            }
        }
        return true;
    }
}

