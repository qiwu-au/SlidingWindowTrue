import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

/**********************
Written by Qi Wu, Department of Biomedicine, Aarhus University
Date: 2020-09-11
**********************/

public class SlidingWindowTrueWriter {
private BufferedReader br;

/**********************
Definition
file: file directory where the input file is placed
THe following code was tested on "KTEA_proteome.txt" as input file, downloaded from https://esbl.nhlbi.nih.gov/KTEA/
**********************/

	public SlidingWindowTrueWriter(String file) throws IOException{
		this.br = new BufferedReader(new FileReader(file));
	}

/**********************
Definition
output: file directory where the output file is placed
numberOfSegment: number of columns with expression values, in the case it is 14
minimumAbsoluteValue: the minimum expression value that is deemed as valid quantification
foldChangeThres: the fold change threshold used for intra-segment expression comparisons
**********************/
	
	public void write (String output, int numberOfSegment, double minimumAbsoluteValue, double foldChangeThres) throws IOException{
		PrintWriter pw = new PrintWriter(output);
		br.readLine();
		String line = null;
		while((line = br.readLine()) != null){	
			
/**********************
Load the data. 
The following code works for the "KTEA_proteome.txt" downloaded from https://esbl.nhlbi.nih.gov/KTEA/
Possible changes include: tab delimiter (comma or tab), number of tabs before the expression values (in this case it is 3)
**********************/
			
			String[] str = line.split("\t");
			double[] valueArray = new double[numberOfSegment];
			for(int i=0; i<numberOfSegment; i++){
				valueArray[i] = Double.parseDouble(str[i+3].trim());
			}
			
/**********************
Pick out single segment that has a higher expression value 
(above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/
					
			boolean[] singleBooleanArray = new boolean [numberOfSegment];
			for(int i=0; i<numberOfSegment; i++){
				int count1=0;
				if(valueArray[i]>minimumAbsoluteValue){
					for(int j=0; j<numberOfSegment; j++){
						if(j!=i){
							if(valueArray[i]>foldChangeThres*valueArray[j]){
								count1++;
							}
						}
					}
					if(count1 == numberOfSegment-1){
						singleBooleanArray[i] = true;
					}
				}
			}
			
/**********************
Pick out continuous double segments (intra-segment expression difference is <"foldChangeThres" times) 
that has higher expression values (both above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/			
					
			boolean[] doubleBooleanArray = new boolean [numberOfSegment-1];
			for(int i=0; i<numberOfSegment-1; i++){
				int count2=0;
				double[] tempArray = new double[2];
				tempArray[0] = valueArray[i];
				tempArray[1] = valueArray[i+1];
				Arrays.sort(tempArray);
				if(tempArray[0]>minimumAbsoluteValue){
					if(tempArray[1]<foldChangeThres*tempArray[0]){
						for(int j=0; j<numberOfSegment; j++){
							if(j!=i && j!=i+1){
								if(tempArray[0]>foldChangeThres*valueArray[j]){
									count2++;
								}
							}	
						}
					}
				}
				if(count2 == numberOfSegment-2){
					doubleBooleanArray[i] = true;
				}
			}

/**********************
Pick out continuous three segments (intra-segment expression differences are <"foldChangeThres" times) 
that has higher expression values (all above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/	
		
			boolean[] threeBooleanArray = new boolean [numberOfSegment-2];
			for(int i=0; i<numberOfSegment-2; i++){
				int count3=0;
				double[] tempArray = new double[3];
				tempArray[0] = valueArray[i];
				tempArray[1] = valueArray[i+1];
				tempArray[2] = valueArray[i+2];
				Arrays.sort(tempArray);
				if(tempArray[0]>minimumAbsoluteValue){
					if(tempArray[2]<foldChangeThres*tempArray[0]){
						for(int j=0; j<numberOfSegment; j++){
							if(j!=i && j!=i+1 && j!=i+2){
								if(tempArray[0]>foldChangeThres*valueArray[j]){
									count3++;
								}
							}	
						}
					}
				}
				if(count3 == numberOfSegment-3){
					threeBooleanArray[i] = true;
				}
			}
			
/**********************
Pick out continuous four segments (intra-segment expression differences are <"foldChangeThres" times) 
that has higher expression values (all above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/		

			boolean[] fourBooleanArray = new boolean [numberOfSegment-3];
			for(int i=0; i<numberOfSegment-3; i++){
				int count4=0;
				double[] tempArray = new double[4];
				tempArray[0] = valueArray[i];
				tempArray[1] = valueArray[i+1];
				tempArray[2] = valueArray[i+2];
				tempArray[3] = valueArray[i+3];
				Arrays.sort(tempArray);
				if(tempArray[0]>minimumAbsoluteValue){
					if(tempArray[3]<foldChangeThres*tempArray[0]){
						for(int j=0; j<numberOfSegment; j++){
							if(j!=i && j!=i+1 && j!=i+2 && j!=i+3){
								if(tempArray[0]>foldChangeThres*valueArray[j]){
									count4++;
								}
							}	
						}
					}
				}
				if(count4 == numberOfSegment-4){
					fourBooleanArray[i] = true;
				}
			}
			
/**********************
Pick out continuous five segments (intra-segment expression differences are <"foldChangeThres" times) 
that has higher expression values (all above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/	

			boolean[] fiveBooleanArray = new boolean [numberOfSegment-4];
			for(int i=0; i<numberOfSegment-4; i++){
				int count5=0;
				double[] tempArray = new double[5];
				tempArray[0] = valueArray[i];
				tempArray[1] = valueArray[i+1];
				tempArray[2] = valueArray[i+2];
				tempArray[3] = valueArray[i+3];
				tempArray[4] = valueArray[i+4];
				Arrays.sort(tempArray);
				if(tempArray[0]>minimumAbsoluteValue){
					if(tempArray[4]<foldChangeThres*tempArray[0]){
						for(int j=0; j<numberOfSegment; j++){
							if(j!=i && j!=i+1 && j!=i+2 && j!=i+3 && j!=i+4){
								if(tempArray[0]>foldChangeThres*valueArray[j]){
									count5++;
								}
							}	
						}
					}
				}
				if(count5 == numberOfSegment-5){
					fiveBooleanArray[i] = true;
				}
			}
			
/**********************
Pick out continuous six segments (intra-segment expression differences are <"foldChangeThres" times) 
that has higher expression values (all above the "minimumAbsoluteValue" and >"foldChangeThres" times higher) than all other segments. 
**********************/	

			boolean[] sixBooleanArray = new boolean [numberOfSegment-5];
			for(int i=0; i<numberOfSegment-5; i++){
				int count6=0;
				double[] tempArray = new double[6];
				tempArray[0] = valueArray[i];
				tempArray[1] = valueArray[i+1];
				tempArray[2] = valueArray[i+2];
				tempArray[3] = valueArray[i+3];
				tempArray[4] = valueArray[i+4];
				tempArray[5] = valueArray[i+5];
				Arrays.sort(tempArray);
				if(tempArray[0]>minimumAbsoluteValue){
					if(tempArray[5]<foldChangeThres*tempArray[0]){
						for(int j=0; j<numberOfSegment; j++){
							if(j!=i && j!=i+1 && j!=i+2 && j!=i+3 && j!=i+4 && j!=i+5){
								if(tempArray[0]>foldChangeThres*valueArray[j]){
									count6++;
								}
							}	
						}
					}
				}
				if(count6 == numberOfSegment-6){
					sixBooleanArray[i] = true;
				}
			}

/**********************
If there is a need to go above six continuous segments,the code can be easily adjusted based on examples set above.
**********************/				

/**********************
Print out the output file.
The content for the final printout can be adjusted based on user need.
**********************/		
			
			pw.print(line);
			
			for(int i=0; i<singleBooleanArray.length; i++){
				pw.print("\t"+singleBooleanArray[i]);
			for(int i=0; i<doubleBooleanArray.length; i++){
				pw.print("\t"+doubleBooleanArray[i]);
			}
			for(int i=0; i<threeBooleanArray.length; i++){
				pw.print("\t"+threeBooleanArray[i]);
			}
			for(int i=0; i<fourBooleanArray.length; i++){
				pw.print("\t"+fourBooleanArray[i]);
			}			
			for(int i=0; i<fiveBooleanArray.length; i++){
				pw.print("\t"+fiveBooleanArray[i]);
			}			
			for(int i=0; i<sixBooleanArray.length; i++){
				pw.print("\t"+sixBooleanArray[i]);
			}			
			pw.print("\n");
		}				
		pw.close();
	}
	
	public static void main(String[] args) throws IOException{
		SlidingWindowTrueWriter swtw = new SlidingWindowTrueWriter("U:\\KTEA_proteome.txt");
		swtw.write("U:\\KTEA_proteome_boolean.txt",14,10,10);
	}

}
