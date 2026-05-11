import edu.rice.melts.MELTS.*;

public class MELTSengine {
    private MELTS melts;

    public MELTSengine() {
        melts = new MELTS();
    }

    public void setBulkComposition(String oxide, double wt) {
        melts.setBulkComposition(oxide, wt);
    }

    public void setPressure(double Pbar) {
        melts.setPressure(Pbar);
    }

    public void setTemperature(double Tc) {
        melts.setTemperature(Tc);
    }

    public void setOxygenFugacity(String buffer) {
        melts.setOxygenFugacity(buffer);
    }

    public void calcEquilibriumState() {
        melts.calcEquilibriumState();
    }

    public String[] getSolidNames() {
        return melts.getSolidNames();
    }

    public String[] getLiquidNames() {
        return melts.getLiquidNames();
    }

    public double getPhaseMass(String name) {
        return melts.getPhaseMass(name);
    }
}
