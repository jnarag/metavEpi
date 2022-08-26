import java.util.ArrayList;
import java.util.List;

public class infections {


    List<Integer> id = null;
    List<Double> birth = null;
    List<Double> death = null;

    public infections(int size) {
        id = new ArrayList<>(size);
        birth = new ArrayList<>(size);
        death = new ArrayList<>(size);
    }

    public infections() {
        id = new ArrayList<>();
        birth = new ArrayList<>();
        death = new ArrayList<>();
    }

    public List<Integer> getId() {
        return id;
    }

    public void setId(List<Integer> id) {
        this.id = id;
    }

    public double getBirth(int index) {
        return birth.get(index);
    }

    public void setBirth(int index, double t) {
        this.birth.set(index, t);
    }

    public void setBirth(Integer infection, double t) {
        this.birth.set(this.id.indexOf(infection), t);
    }

    public Double getDeath(int index) {

        return death.get(index);
    }

    public void setDeath(int index, double t) {

        this.death.set(index, t);
    }

    public void setDeath(Integer infection, double t) {

        this.death.set(this.id.indexOf(infection), t);
    }


    public void logInfection(int id, double birth, double death) {

        this.id.add(id);
        this.birth.add(birth);
        this.death.add(death);

    }


}
