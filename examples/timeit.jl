macro timeit(ex,name)
    t, i = gensym(2)
    quote
        $t = 0
        for $i=1:20
            $t += @elapsed $ex
        end
        $t = $t/20
        printf("%s took %.1f ms (avg over 20 runs)\n", $name, $t*1000)
    end
end
