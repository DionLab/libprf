package PfTools;

public class Matrix {
	static {
		try {
			//String OS = System.getProperty("os.name").toLowerCase();
			System.loadLibrary("prf");
		} catch (Exception e) {
			System.exit(-1);
		}
	}
  
  public Matrix( Profile prf, String Sequence) {
		SequenceLength = Sequence.length();
		matrixPointer = 0;
		profilePointer = prf.profilePointer;
		workPointer = prf.workPointer;
		if (! buildMatrix(Sequence)) {
			throw new RuntimeException("Error computing matrix");
		}
  }
  
	@Override
	public void finalize() {
		freeMatrix();
	}
  
	private long matrixPointer;
	private long SequenceLength;
	private long profilePointer;
	private long workPointer;
	
	native private void freeMatrix();
	native private boolean buildMatrix(String sequence);
	native public void dumpTabulated(String FileName);
}
