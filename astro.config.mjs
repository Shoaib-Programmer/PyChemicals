// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";

export default defineConfig({
  image: {
    service: { entrypoint: "astro/assets/services/noop" },
  },
  site: "https://Shoaib-Programmer.github.io/PyChemicals",
  base: "/PyChemicals/",
  integrations: [
    starlight({
      title: "PyChemicals",
      social: {
        github: "https://github.com/Shoaib-Programmer/PyChemicals",
      },
      sidebar: [
        {
          label: "Guides",
          items: [{ label: "Example Guide", link: "/guides/example/" }],
        },
        {
          label: "Reference",
          autogenerate: { directory: "reference" },
        },
      ],
      customCss: ["./src/styles/globals.css"],
      defaultLocale: "en",
      locales: {
        en: { label: "English" },
        es: { label: "Español" },
        fr: { label: "Français" },
      },
    }),
  ],
});
