#use "topfind";;
#thread;;
#require "ketrew";;
#require "coclobas.ketrew_backend";;
open Printf;;

let () = Unix.putenv "KETREW_CONFIGURATION" (Sys.getenv "KETREW_ROOT" ^ "/configuration.ml");;

let run_2_commands_with_python_method cmd1 =
  let open Ketrew.EDSL in
  let node1 =
    workflow_node without_product
      ~name:"Pyd:Prep"
      ~make:(
        Coclobas_ketrew_backend.Plugin.local_docker_program
          ~tmp_dir:"/tmp/secotrec-local-shared-temp"
          ~base_url:"http://coclo:8082"
          ~image:"py3.5-v1"
          ~volume_mounts:[
            `Local ("/Users/eczech/projects/hammer/repos/cartpipeline/", "/tmp/cartpipeline")
          ]
          Program.(
            sh "export PYTHONPATH=\"/tmp/cartpipeline/pyhpa:/tmp/cartpipeline/pycgds:$PYTHONPATH\"" &&
            sh "python /tmp/cartpipeline/pyhpa/script/gene_selector.py --output /tmp/cartpipeline/test.csv"
          )
          (*(pypgrm cmd1)*)
      )
  in
  Ketrew.Client.submit_workflow node1;;

let () = run_2_commands_with_python_method "ls";;
