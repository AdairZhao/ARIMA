export type Itinerary = {
  slug: string;
  title: string;
  days: number;
  summary: string;
  image: string;
  imageAlt: string;
  route: string[];
  daysPlan: { day: string; focus: string; notes: string }[];
};

export const itineraries: Itinerary[] = [
  {
    slug: "classic-10",
    title: "Classic first trip · 10 days",
    days: 10,
    summary:
      "Beijing → Xi'an → Shanghai. History, Wall, Warriors, then a soft landing in China’s most international city.",
    image: "/images/palace.jpg",
    imageAlt: "The Forbidden City in Beijing",
    route: ["Beijing", "Xi'an", "Shanghai"],
    daysPlan: [
      {
        day: "Days 1–4",
        focus: "Beijing",
        notes:
          "Forbidden City, Temple of Heaven, hutongs, Mutianyu Wall. Jet-lag buffer on day 1.",
      },
      {
        day: "Days 5–6",
        focus: "Xi'an",
        notes:
          "Evening arrival by high-speed rail. Warriors + wall + Muslim Quarter food.",
      },
      {
        day: "Days 7–10",
        focus: "Shanghai",
        notes:
          "Bund, Concession walks, dumplings, optional water town half-day. Fly out from PVG/SHA.",
      },
    ],
  },
];

export function getItinerary(slug: string) {
  return itineraries.find((i) => i.slug === slug);
}
