function spanLift = findSpanLift(results,geo,state)

lemma=size(results.F);
for i=1:lemma(1)
    B2WTransform=[cos(state.betha)*cos(state.alpha(i)),        -sin(state.betha),          cos(state.betha)*sin(state.alpha(i)) ;...
              cos(state.alpha(i))*sin(state.betha),         cos(state.betha),          sin(state.betha)*sin(state.alpha(i)) ;...
                              -sin(state.alpha(i)),                        0,                           cos(state.alpha(i))];
    lemma4(i,:)=B2WTransform*results.F(i,:)';                                         
	forceLift(i)=lemma4(i,3);                     %Lift on each panel
end

spanLift = zeros(geo.Wings(1).wing.SegNum,1);
for i = 1:geo.Wings(1).wing.SegNum
    spanLift(i,1) = sum(forceLift(i:geo.Wings(1).wing.SegNum:geo.Wings(1).wing.cordNum*geo.Wings(1).wing.SegNum));
end