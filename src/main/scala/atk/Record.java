package atk;

/**
 * Copied from Thomas Abeel's ATK on 29/11/2017
 */

public class Record {
    private String sequence;

    private String description;

    public Record(String description, String sequence) {
        this.sequence = sequence;
        this.description = description;
    }

    public boolean equals(Object o) {
        Record x = (Record) o;
        return x.getDescription().equals(this.getDescription()) && x.getSequence().equals(this.getSequence());

    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

}