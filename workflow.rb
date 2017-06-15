
require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/util/R'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/TCPA'

TCPA_BATCHES = TCPA.cell_lines.produce.glob("*").collect{|f| File.basename(f) }
module TCPA
  extend Workflow

  ASMUND_ANTIBODIES = TCPA.Asmund.tsv.fields
  TCPA_BATCHES = TCPA.cell_lines.produce.glob("*").collect{|f| File.basename(f) }
  TCPA_ANTIBODIES = TCPA.cell_lines[TCPA_BATCHES.first].tsv.fields
  ANTIBODIES = TCPA_ANTIBODIES & ASMUND_ANTIBODIES

  helper :rppa_matrix do
    TCPA.all_samples_fixed.find
  end

  input :batch, :select, "Select batch", nil, :select_options => TCPA_BATCHES
  task :normalized_tcpa => :tsv do |batch|
    tsv = TCPA.cell_lines[batch].tsv
    new = tsv.R <<-EOF
data.new = data.frame(t(apply(data,1,function(x){sort(x, index.return=T)$ix / length(x)})))
rownames(data.new) = rownames(data)
names(data.new) = names(data)
data = data.new
    EOF
    new.key_field = "Sample"
    new.cast = :to_f
    new
  end

  task :normalized_asmund => :tsv do 
    tsv = TCPA.Asmund.tsv
    new = tsv.R <<-EOF
data.new = data.frame(t(apply(data,1,function(x){sort(x, index.return=T)$ix / length(x)})))
rownames(data.new) = rownames(data)
names(data.new) = names(data)
data = data.new
    EOF
    new.key_field = "Sample"
    new.cast = :to_f
    new
  end

  dep :normalized_tcpa
  dep :normalized_asmund
  task :antibody_analysis => :tsv do 
    asmund = step(:normalized_asmund).load
    tcpa = step(:normalized_tcpa).load
    controls = asmund.keys.select{|k| k =~ /Ctr/}

    tsv = TSV.setup({}, :key_field => "Antibody", :fields => ["Asmund mean", "TCPA mean", "High", "P-value"], :type => :list)
    ANTIBODIES.each do |antibody|
      begin
        asmund_values = asmund.column(antibody)
        tcpa_values = tcpa.column(antibody)
        tcpa_scores = tcpa_values.values
        asmund_scores = asmund_values.values_at(*controls)
        cmd = "t.test(#{R.ruby2R asmund_scores},#{R.ruby2R tcpa_scores})$p.value"
        p_value = R.eval cmd
        mean_tcpa = Misc.mean(tcpa_scores)
        mean_asmund = Misc.mean(asmund_scores)
        tsv[antibody] = [mean_asmund, mean_tcpa, (mean_asmund > mean_tcpa).to_s, p_value]
      rescue
        Log.exception $!
        next
      end
    end

    tsv
  end

  task :ranks => :tsv do 
    rppa_matrix.tsv.R <<-EOF
d = apply(data,1,function(x){ rbbt.ranks(x) })
d = data.frame(t(d))
rownames(d) <- rownames(data)
names(d) <- names(data)
data = d
    EOF
  end

  dep :ranks
  input :sample, :select, "Sample description", nil, :select_options => TCPA.all_samples_fixed.tsv.keys.sort
  input :fold_change, :select, "Standard deviations from mean threshold", 1, :select_options => [0,0.1,0.5,1,1.5,2]
  task :sample_activation => :tsv do |sample,fold_change|
    fold_change = fold_change.to_f
    ranks = step(:ranks).load
    other_samples = ranks.keys - [sample]
    sample_ranks = ranks[sample]
    ab_level = TSV.setup({}, :key_field => "Antibody", :fields => ["Abundance", "Value", "Mean", "SD"], :type => :list)
    ranks.fields.each do |antibody|
      column = ranks.column(antibody)
      value = column[sample].to_f
      other = column.values.compact.reject{|v| v.empty? or v == "NA"}.collect{|v| v.to_f}
      mean = Misc.mean(other)
      sd = Misc.sd(other)
      active = case 
               when value > mean + fold_change * sd
                 "1"
               when value < mean - fold_change * sd
                 "0"
               else
                 "-"
               end
      active = nil if value == 0
      ab_level[antibody] = [active, value, mean,sd]
    end
    ab_level
  end

  input :fold_change, :select, "Standard deviations from mean threshold", 1, :select_options => [0,0.1,0.5,1,1.5,2]
  task :activations => :tsv do |fold_change|
    samples = rppa_matrix.tsv.keys
    result = nil
    TSV.traverse samples, :type => :array, :bar => self.progress_bar("Processing samples") do |sample|
      tsv = TCPA.job(:sample_activation, nil, :sample => sample, :fold_change => fold_change).run
      tsv = tsv.column("Abundance")
      tsv.fields = [sample]

      if result.nil?
        result = tsv
      else
        result = result.attach tsv, :complete => true
      end
    end
    result.to_s
  end


  export_asynchronous :ranks, :sample_activation, :activations
end

#require 'rbbt/knowledge_base/TCPA'
#require 'TCPA/tasks/activity.rb'
#require 'rbbt/entity/TCPA'

