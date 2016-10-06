#KnowledgeBaseRESTHelpers.add_syndication :TCPA, TCPA.knowledge_base


module Object::Antibody
  extend Entity
  include Entity::REST

end

get '/antibody' do
  template_render('antibody', @clean_params, "Combination", :cache_type => :sync)
end

get '/correlation' do
  template_render('correlation', @clean_params, "Combination", :cache_type => :sync)
end
