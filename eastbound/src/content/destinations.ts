export type Destination = {
  slug: string;
  name: string;
  region: string;
  blurb: string;
  bestFor: string;
  days: string;
  image: string;
  imageAlt: string;
  see: string[];
  eat: string[];
  move: string[];
  stay: string;
};

export const destinations: Destination[] = [
  {
    slug: "beijing",
    name: "Beijing",
    region: "North",
    blurb:
      "Imperial courtyards, hutong alleys, and the Great Wall within a day trip — China’s political and cultural capital.",
    bestFor: "History first-timers",
    days: "4–5 days",
    image: "/images/great-wall.jpg",
    imageAlt: "The Great Wall of China along misty ridgelines",
    see: [
      "Forbidden City (book passport tickets ahead)",
      "Mutianyu Great Wall (less crowded than Badaling)",
      "Temple of Heaven at opening hour",
      "Hutong walk around Nanluoguxiang — then leave the main street",
    ],
    eat: [
      "Peking duck (look for local lunch sets)",
      "Zhajiangmian in a hutong shop",
      "Night snacks at Dongzhimen or Huguo Temple area",
    ],
    move: [
      "Metro covers most sights; get a transit QR via Alipay",
      "DiDi for late nights and Wall transfers",
      "Allow buffer time — security checks are common",
    ],
    stay: "Dongcheng or Chaoyang for metro access to major sights.",
  },
  {
    slug: "shanghai",
    name: "Shanghai",
    region: "East",
    blurb:
      "Skyline energy, concession-era streets, and easy day trips — the softest landing for many first visits.",
    bestFor: "City + food",
    days: "3–4 days",
    image: "/images/shanghai.jpg",
    imageAlt: "Shanghai Pudong skyline across the river",
    see: [
      "The Bund at blue hour",
      "Yu Garden + Old Street (go early)",
      "French Concession tree-lined walks",
      "Shanghai Tower observation (clear-day tickets)",
    ],
    eat: [
      "Xiaolongbao (soup dumplings) — share orders",
      "Shengjianbao breakfast",
      "Riverfront dinner with a reservation",
    ],
    move: [
      "Metro is excellent; Maglev optional from Pudong airport",
      "Walk the Concession — taxis stall in peak traffic",
    ],
    stay: "Huangpu or Jing’an for walkability; Pudong if you want towers.",
  },
  {
    slug: "xian",
    name: "Xi'an",
    region: "Northwest",
    blurb:
      "Terracotta Warriors, city walls, and Muslim Quarter spices — the Silk Road still starts here.",
    bestFor: "History + night food",
    days: "2–3 days",
    image: "/images/xian.jpg",
    imageAlt: "Rows of Terracotta Warriors in Xi'an",
    see: [
      "Terracotta Army (book shuttle or organized transfer)",
      "Bike or walk the Ming city wall",
      "Giant Wild Goose Pagoda area at night",
    ],
    eat: [
      "Roujiamo (Chinese “burger”)",
      "Yangrou paomo — ask for tourist-friendly broth tips",
      "Muslim Quarter skewers after dark",
    ],
    move: [
      "Metro to wall & pagoda; Warriors need a coach or DiDi",
      "High-speed rail from Beijing or Chengdu is straightforward",
    ],
    stay: "Inside or near the city wall for evening food walks.",
  },
  {
    slug: "chengdu",
    name: "Chengdu",
    region: "Southwest",
    blurb:
      "Pandas, teahouses, and Sichuan heat — a slower rhythm that still feeds you well.",
    bestFor: "Food + pandas",
    days: "3–4 days",
    image: "/images/chengdu.jpg",
    imageAlt: "Giant pandas at a Chengdu research center",
    see: [
      "Chengdu Research Base of Giant Panda Breeding (early entry)",
      "People’s Park teahouse people-watching",
      "Jinli / Kuanzhai Alley — brief visits, then find quieter streets",
      "Optional Leshan or Mount Qingcheng day trip",
    ],
    eat: [
      "Hotpot — ask for mild / foreigner-friendly spice levels",
      "Dan dan noodles and mapo tofu at local shops",
      "Late-night skewers",
    ],
    move: [
      "Metro is growing fast; DiDi fills gaps",
      "Panda base: go at opening, then nap — Chengdu afternoons are for tea",
    ],
    stay: "Jinjiang or Qingyang for central metro and food streets.",
  },
];

export function getDestination(slug: string) {
  return destinations.find((d) => d.slug === slug);
}
