import org.apache.groovy.json.internal.IO;

import java.io.*;
import java.util.stream.IntStream;
import java.util.Random;

public class VCFBuilder {
    private String header = "##fileformat=VCFv4.2\n" +
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
            "##samtoolsVersion=1.12+htslib-1.12\n" +
            "##samtoolsCommand=samtools mpileup -d 250 -ugf Homo_sapiens.GRCh37.dna.primary_assembly.gz read.sorted.bam\n" +
            "##reference=file:Homo_sapiens.GRCh37.dna.primary_assembly.gz\n" +
            "##contig=<ID=1,length=249250621>\n" +
            "##contig=<ID=10,length=135534747>\n" +
            "##contig=<ID=11,length=135006516>\n" +
            "##contig=<ID=12,length=133851895>\n" +
            "##contig=<ID=13,length=115169878>\n" +
            "##contig=<ID=14,length=107349540>\n" +
            "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n" +
            "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n" +
            "##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">\n" +
            "##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">\n" +
            "##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">\n" +
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            "##bcftools_callVersion=1.12+htslib-1.12\n" +
            "##bcftools_callCommand=call -vmO z -o 'variants.vcf.gz' 'read.pileup'; Date=Sun May 27 15:04:45 2021\n";
    private String rathausBar = "#CHROM\tPOS\t\t\t\tREF\t    ALT";
    private String hospitalBar = "##CHROM\tPOS\t\t\t\tID\t\tREF\t    ALT\t\tQUAL\t\tFILTER\tINFO\t\tFORMAT";
    private String[] qual = {"23.3425", "24.2332", "53.3321", "54.2354", "43.2355", "41.2345", "16.2155", "32.1252", "21.2245"};
    private String[] info = {
            "DP=16;VDB=0.907611;SGB=-0.636426;RPB=0.966012;MQB=0.428703;BQB=0.0628765;MQ0F=0;AC=1;AN=2;DP4=7,0,7,0;MQ=49",
            "DP=17;VDB=0.539658;SGB=-0.636426;RPB=0.887766;MQB=0.621145;BQB=0.708895;MQ0F=0;AC=1;AN=2;DP4=10,0,7,0;MQ=54",
            "DP=16;VDB=0.0485232;SGB=-0.616816;RPB=0.686279;MQB=0.863243;BQB=0.0253122;MQ0F=0;AC=1;AN=2;DP4=10,0,6,0;MQ=57",
            "DP=16;VDB=0.0485232;SGB=-0.616816;RPB=0.686279;MQB=0.863243;BQB=0.0292791;MQ0F=0;AC=1;AN=2;DP4=10,0,6,0;MQ=57",
            "DP=28;VDB=0.877004;SGB=-0.680642;RPB=0.877755;MQB=0.933359;BQB=0.0384;MQ0F=0;AC=1;AN=2;DP4=16,0,12,0;MQ=58",
            "DP=29;VDB=0.67865;SGB=-0.676189;RPB=0.923174;MQB=1;BQB=0.628158;MQ0F=0;AC=1;AN=2;DP4=17,0,11,0;MQ=60",
            "DP=49;VDB=0.245012;SGB=-0.692976;RPB=0.976675;MQB=3.31401e-07;BQB=1.09401e-05;MQ0F=0.0204082;AC=1;AN=2;DP4=22,0,26,0;MQ=37",
            "DP=26;VDB=0.120141;SGB=-0.692976;MQ0F=0.0769231;AC=2;AN=2;DP4=0,0,26,0;MQ=12",
            "DP=5;VDB=0.309755;SGB=-0.511536;RPB=0.333333;MQB=1;BQB=0;MQ0F=0;AC=1;AN=2;DP4=2,0,3,0;MQ=60"
    };
    private char[] base = {'A','G','C','T'};
    private int chrom = 3;

    public void writeVCFHospital(String path) {
        File file = new File(path);
        String output = header+hospitalBar;
        for(int i = 10035634; i <= 10035734; i++) {
            int randBase = new Random().nextInt(base.length);
            int randRef = new Random().nextInt(base.length);
            int randQual = new Random().nextInt(qual.length);
            int randInfo = new Random().nextInt(info.length);
            output = output+"\n"+chrom+"\t    "+i+"\t    .\t    "+base[randBase]+"\t    "+base[randRef]+"\t    "+qual[randQual]+"\t    .\t    "+info[randInfo]+"\t    GT:PL";
        }
        try {
            FileWriter writer = new FileWriter(file);
            writer.write(output);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void writeVCF(String path) {
        File file = new File(path);
        try {
            file.createNewFile();
        } catch (IOException io) {
            io.printStackTrace();
        }
        String output = header+rathausBar;
        for(int i = 10035680; i <= 10035700; i++) {
            int randBase = new Random().nextInt(base.length);
            int randRef = new Random().nextInt(base.length);
            output = output+"\n"+chrom+"\t    "+i+"\t    "+base[randBase]+"\t    "+base[randRef];
        }
        try {
            FileWriter writer = new FileWriter(file);
            writer.write(output);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        VCFBuilder builder = new VCFBuilder();
        //builder.writeVCFHospital(".\\VCFs\\Patient01");
        for(int i = 1; i < 11; i++) {
            builder.writeVCF(".\\VCFs\\Citizen"+i+".csv");
        }
        for(int i = 1; i < 11; i++) {
            builder.writeVCFHospital(".\\VCFs\\Patient"+i+".csv");
        }
    }
}
