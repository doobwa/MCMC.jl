macro timeit(ex,name)
    t, i = gensym(2)
    quote
        $t =  @elapsed $ex
        printf("%s took %.1f ms\n", $name, $t*1000)
    end
end
