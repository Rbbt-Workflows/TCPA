require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module TCPA
  extend Resource
  self.subdir = 'share/databases/TCPA'

  #self.search_paths = {}
  #self.search_paths[:default] = :workflow
  #self.search_paths[:workflow] = Path.caller_lib_dir 

  def self.organism
    Organism.default_code "Hsa"
  end

  TCPA.claim TCPA.MDACC["data.zip"], :proc do |filename|
    raise "Download the MDACC (all of it) from http://app1.bioinformatics.mdanderson.org/tcpa/_design/basic/download.html and place it in #{filename}"
  end


  TCPA.claim TCPA.antibodies, :proc do
    tsv = Rbbt.data["antibody_info"].find(:lib).tsv :type => :list, :indentifiers => false
    tsv.add_field "Feature" do |name, values|
      code = values.first
      if code =~ /cleav/
        "cleaved"
      else
        code = code.sub(/-([VCRM-]+)_(GBL.*)/,'_\2-\1')
        code = code.sub(/(\d)_([pYTS]\d+)/,'\1-\2')
        _g, site, abc = code.split "_"
        site || "abundance"
      end
    end

    tsv.to_s
  end

  TCPA.claim TCPA.identifiers, :proc do |filename|
    zipfile = TCPA.MDACC["data.zip"].produce.find
    source_dir = Rbbt.tmp[".source"].find
    if not File.exists? source_dir
      FileUtils.mkdir_p source_dir
      CMD.cmd("cd '#{source_dir}' && unzip '#{zipfile}'")
    end

    require 'rbbt/ner/rnorm'
    Log.info "Prepareing norm"
    ids = Rbbt.data["antibody_identifiers"].tsv
    ids2 = Rbbt.data["antibody_identifiers_final"].tsv
    ids_tmp = ids.change_key("Ab_Full name", :identifiers => ids2).to_flat
    
    ids = {}
    ids_tmp.each do |k,values|
      ids[k.gsub(/ +/,'')] = values.collect{|v| v.gsub(/ +/,'') }
    end
    ids_tmp.annotate(ids)
    @norm = Normalizer.new ids
    Log.info "Prepared norm"
    final = TSV.setup(ids.keys, :key_field => ab_info.key_field, :fields => [], :type => :double)
    source_dir.glob("*.csv").each do |csvfile|
      Log.info "Processing: #{ csvfile }"
      tsv = csvfile.tsv(:sep => ',', :header_hash => '', :fix => Proc.new{|line| line.gsub(/"|\.txt/,'')}, :type => :list, :key_field => '"Sample description"', :cast => :to_f)

      tsv.fields = tsv.fields.collect{|f| f.gsub('"','').gsub('.txt','').gsub(/ +/,'')}
      tsv.key_field = tsv.key_field.gsub('"','')

      good_fields = tsv.fields[7..-1]
      tsv = tsv.slice(good_fields)
      trans = TSV.setup({}, :key_field => ab_info.key_field, :fields => [File.basename(csvfile).sub(/.csv/,'')], :type => :flat)
      missing = []
      tsv.fields.each do |field|
        ffield = field.sub(/(-[A-Z])*_G..\d+$/,'')
        match = @norm.resolve(ffield).first
        if match
          trans[match] ||= []
          trans[match] << field
        else
          missing << field
        end
      end

      final = final.attach trans
    end

    Open.write(filename, final.to_s)
    nil
  end

  TCPA.claim TCPA.NCI60, :proc do |filename|

    translations = Rbbt.data["abs_S29_L3_S29_names.tsv"].index :target => "Ab Full name_batch_29"
    zipfile = TCPA.MDACC["data.zip"].find
    tsv = nil
    TmpFile.with_file do |source_dir|
      Path.setup(source_dir)
      if not File.exists? source_dir
        FileUtils.mkdir_p source_dir
        CMD.cmd("cd '#{source_dir}' && unzip '#{zipfile}'")
      end
      csvfile = source_dir.glob("*.csv").select{|f| f =~ /S29/ and not f =~ /TREAT/  }.first

      tsv = csvfile.tsv(:sep => ',', :header_hash => '', :fix => Proc.new{|line| line.gsub('"','')}, :type => :double, :merge => true, :key_field => '"Sample description"', :cast => :to_f)
    end

    tsv.fields = tsv.fields.collect{|f| f.gsub('"','') }

    good_fields = tsv.fields[7..-1]
    tsv = tsv.slice(good_fields)
    new_fields = translations.values_at(*tsv.fields)
    tsv.fields = new_fields
    tsv.key_field = tsv.key_field.gsub('"','')

    tsv.add_field "Replicate" do |sample,values|
      (1..values.first.length).to_a
    end

    new_tsv = tsv.unzip("Replicate", true)

    Open.write(filename, new_tsv.to_s)
    nil
  end

  TCPA.claim TCPA.cell_lines, :proc do |directory|
    zipfile = TCPA.MDACC["data.zip"].produce.find
    source_dir = directory[".source"].find
    FileUtils.mkdir_p source_dir
    CMD.cmd("cd '#{source_dir}' && unzip '#{zipfile}'")
    source_dir.glob("*.csv").each do |csvfile|
      tsv = csvfile.tsv(:sep => ',', :header_hash => '', :fix => Proc.new{|line| line.gsub(/"|\.txt/,'')}, :type => :list, :key_field => '"Sample description"', :cast => :to_f)

      tsv.fields = tsv.fields.collect{|f| f.gsub('"','').gsub('.txt','')}
      tsv.key_field = tsv.key_field.gsub('"','')

      good_fields = tsv.fields[7..-1]
      tsv = tsv.slice(good_fields)

      #tsv.fields = tsv.fields.collect{|f| f.sub(/(-[MVRQCENA]+)+[-_]+GB[L,J]\d+/,'').sub(/-(\d+)/,'\1').sub('_','-').upcase.sub('B CATE', 'BETA CATE')}

      Open.write(directory[File.basename(csvfile).sub(/\.csv$/,'')], tsv.to_s)
    end
  end

  TCPA.claim TCPA.all_samples, :proc do |filename|
    all = nil
    name_index = TCPA.identifiers.index
    all_abs = name_index.values.flatten.compact.uniq.sort
    missing = nil
    #all = TSV.setup({}, :key_field => "Sample description", :fields => , :type => :list, :cast => :to_f)
    TCPA.cell_lines.glob("*").each do |file|
      tsv =TSV.open(file)
      good_fields = []
      fixed_fields = []
      tsv.fields.each do |field|
        next unless name_index.include? field
        good_fields << field
        fixed_fields << name_index[field]
      end
      tsv = tsv.slice good_fields
      tsv.fields = name_index.values_at *good_fields
      missing = all_abs - tsv.fields
      missing.each do |ab|
        tsv.add_field ab do
          nil
        end
      end
      tsv  = tsv.reorder(:key, all_abs)
      if all.nil?
        all = tsv
      else
        tsv.each do |k,v|
          all[k] = v
        end
      end
    end
    all.to_s
  end

  TCPA.claim TCPA.extended_samples, :proc do |filename|
    all = TCPA.all_samples.tsv
    name_index = TCPA.identifiers.index
    all_abs = name_index.values.flatten.compact.uniq.sort
    missing = nil

    tsv = Rbbt.data["Asmund_RPPA_log2norm"].tsv(:type => :list).transpose("Sample description")
    good_fields = []
    fixed_fields = []
    tsv.fields.each do |field|
      next unless name_index.include? field
      good_fields << field
      fixed_fields << name_index[field]
    end
    tsv = tsv.slice good_fields
    tsv.fields = name_index.values_at *good_fields
    missing = all_abs - tsv.fields
    missing.each do |ab|
      tsv.add_field ab do
        nil
      end
    end
    tsv  = tsv.reorder(:key, all_abs)
    tsv.each do |k,v|
      all[k] = v
    end
    all.to_s
  end

  TCPA.claim TCPA.all_samples_fixed, :proc do |filename|
    tsv = nil
    Rbbt.data.latest.glob("*.txt").each do |file|
      new = file.tsv :type => :list, :header_hash => "", :cast => :to_f, :unnamed => true
      if tsv.nil?
        tsv = new
      else

        all_fields = tsv.fields + new.fields
        good_fields = all_fields.uniq.sort
        missing_fields = good_fields - tsv.fields
        missing_fields.each do |field|
          tsv = tsv.add_field(field){ nil }
        end

        missing_fields = good_fields - new.fields
        missing_fields.each do |field|
          new = new.add_field(field){ nil }
        end
        tsv = tsv.reorder(:key, good_fields)
        new = new.reorder(:key, good_fields)
        new.each do |k,v|
          tsv[k] = v
        end
      end
      Log.tsv tsv
    end
    tsv.to_s
  end

  TCPA.claim TCPA.Asmund, :proc do |filename|
    tsv = Rbbt.data["Asmund_RPPA"].tsv(:type => :list, :fix => Proc.new{|l| l.gsub(',','.')}).reorder(5, (9..225).to_a)
    #tsv.fields = tsv.fields.collect{|f| f.sub(/(-[MVRQCENA]+)+[-_]+GB[L,J]\d+/,'').sub(/-(\d+)/,'\1').sub('_','-').upcase.sub('B CATE', 'BETA CATE')}
    tsv.cast = :to_f
    tsv.to_s
  end

  #SAMPLE_FIELDS = ['Sample Source','Sample', 'Order', 'Category', 'Sample description', 'Var.7', 'Revised Antibody Name']

  #TCPA.claim TCPA.RPPA.cell_line_rppa, :proc do 
  #  global_tsv = nil

  #  TCPA.cell_lines.glob("*.csv").each do |file|

  #    # There is a reason to use only this one
  #    next unless file.include? "MDACC-CELLLINE_S29-L3-S29.csv"

  #    tsv = file.tsv(:sep => ',', :header_hash => '', :fix => Proc.new{|line| line.gsub(/"|\.txt/,'')}, :type => :list, :key_field => "Sample")
  #    tsv.fields = tsv.fields.collect{|f| f.gsub('"','').gsub('.txt','')}
  #    tsv.key_field = tsv.key_field.gsub('"','')
  #    global_tsv = global_tsv.nil? ? tsv : global_tsv.merge(tsv)
  #  end

  #  global_tsv.to_s
  #end

  #TCPA.claim TCPA.RPPA.labels, :proc do 
  #  TCPA.RPPA.cell_line_rppa.tsv(:fields => SAMPLE_FIELDS).to_s
  #end

  #TCPA.claim TCPA.RPPA.data, :proc do 
  #  fields = TSV.parse_header(TCPA.RPPA.cell_line_rppa).fields
  #  TCPA.RPPA.cell_line_rppa.tsv(:fields => fields - SAMPLE_FIELDS).transpose("Antibody code").to_s
  #end

  #TCPA.claim TCPA.RPPA.identifiers, :proc do
  #  info = TCPA['.source'].antibody_genes.tsv :key_field => "TCPA antibody name", :fields => ["HGNC"]
  #  info.key_field = "Antibody code"
  #  info.fields = ["Associated Gene Name"]

  #  info.add_field "Clean name" do |code, values|
  #    name = values.first.first
  #    code = code.sub(/-([VCRM-]+)_(GBL.*)/,'_\2-\1')
  #    code = code.sub(/(\d)_([pYTS]\d+)/,'\1-\2')
  #    _g, site, abc = code.split "_"
  #    [name, site, abc].compact * "_"
  #  end

  #  ensembl = Organism.identifiers(TCPA.organism).index :target => "Ensembl Gene ID", :order => true, :persist => true

  #  info.add_field "Ensembl Gene ID" do |code, values|
  #    name = values.first
  #    ensembl.values_at(*name).flatten.uniq
  #  end

  #  info.add_field "PTM" do |code, values|
  #    clean = values["Clean name"].first
  #    g,s,ab = clean.split("_")
  #    ab, s = s, ab if ab.nil?
  #    s || "abundance"
  #  end

  #  info.add_field "Cannonical" do |code, values|
  #    gene = values["Associated Gene Name"].first
  #    ptm = values["PTM"].first
  #    if gene and not gene.empty?
  #      [[gene, ptm] * ":"]
  #    else
  #      []
  #    end
  #  end

  #  info.add_field "Slide code" do |code, values|
  #    clean = values["Clean name"].first
  #    clean.split("_").last
  #  end

  #  info
  #end
end

TCPA.MDACC["data.zip"].produce if __FILE__ == $0
#TCPA.cell_lines.produce if __FILE__ == $0
#TCPA.antibodies.produce(true) if __FILE__ == $0
#TCPA.Asmund.produce(true) if __FILE__ == $0
#TCPA.NCI60.produce(true) if __FILE__ == $0
#TCPA.identifiers.produce(true) if __FILE__ == $0
#TCPA.all_samples.produce(true) if __FILE__ == $0
#TCPA.extended_samples.produce(true) if __FILE__ == $0
#TCPA.all_samples_fixed.produce(true) if __FILE__ == $0
