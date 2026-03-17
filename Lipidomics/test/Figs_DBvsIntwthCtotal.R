
df<-df_with_unsat|>
  filter(!lipid_class1=="FA")|>
  filter(!lipid_class1=="Cer")



ggplot(df, aes(x = DB_total,
                         fill = factor(lipid_class1))) +
  geom_bar(position = position_dodge()) +
  theme_bw(base_size = 12) +
  facet_wrap(~ sname) +
  labs(fill = "DB_total")


ggplot(df, aes(x = DB_total,
               fill = factor(lipid_class1))) +
  geom_bar() +
  theme_bw(base_size = 12) +
  facet_wrap(~ sname) +
  labs(fill = "lipid_class1")

ggplot(df, aes(x = C_total,
               fill = factor(lipid_class1))) +
  geom_bar() +
  theme_bw(base_size = 12) +
  facet_wrap(~ sname) +
  labs(fill = "lipid_class1")


ggplot(df, aes( x= C_total, y = DB_total , fill= lipid_class1,))+
  geom_col(position = position_dodge())+
  theme_bw(base_size = 12)+
  facet_wrap(~sname)


ggplot(df, aes( x= lipid_class1, y = DB_total , fill= sname))+
  geom_col(position = position_dodge())+
  theme_bw(base_size = 12)
  #facet_wrap(~sname)


ggplot(df, aes(x=C_total, y= DB_total, color =lipid_class1))+
  geom_point(alpha=0.6, size = 2)+
  theme_bw(base_size = 12)+
  facet_wrap(~sname)




df_count <- df %>%
  count(sname, lipid_class1, C_total, DB_total, name = "n_lipids")


ggplot(df_count,
       aes(x = C_total,
           y = DB_total,
           color = lipid_class1,
           size = n_lipids)) +
  geom_point(alpha = 0.7) +
  theme_bw(base_size = 12) +
  facet_wrap(~ sname+lipid_class1) +
  scale_size_continuous(range= c(2,6),name = "Number of lipids")


df_count <- df_count |>
  filter(lipid_class1 %in% c("PC", "PE"))
ggplot(df_count,
       aes(x = C_total,
           y = DB_total,
           color = factor(n_lipids),
           size = n_lipids)) +
  geom_point(alpha = 0.7) +
  theme_bw(base_size = 12) +
  facet_wrap(~ sname+lipid_class1, nrow =2) +
  scale_size_continuous(range= c(2,4),name = "Number of lipids")

#==========================
plot<-joined_df|>
  filter(!lipid_class1_int=="FA")|>
  filter(!lipid_class1_int=="Cer")


plot <- plot |>
  filter(lipid_class1_int %in% c("DG", "TG"))

plot <- plot |>
  filter(DB_total >= 6)

ggplot(plot, aes(x=DB_total, y= log2(Intensity),color = factor(C_total)))+
  geom_point(alpha=0.7, size = 2)+
  theme_bw(base_size = 12)+
  facet_wrap(~Sampletype+lipid_class1_int)

ggplot(plot, aes(x=DB_total, y= log2(Intensity),color = factor(C_total)))+
  geom_point(alpha=0.7, size = 2)+
  theme_bw(base_size = 12)+
  facet_wrap(~Sampletype+lipid_class1_int)



