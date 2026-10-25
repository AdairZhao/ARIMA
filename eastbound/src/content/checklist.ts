export type ChecklistItem = {
  id: string;
  title: string;
  detail: string;
  timing: string;
};

export const planChecklist: ChecklistItem[] = [
  {
    id: "visa",
    title: "Visa or transit exemption",
    detail:
      "Confirm whether you need a visa, or if a 144/240-hour transit exemption applies to your passport and routing. Print the latest embassy rules — screenshots alone are not enough at the border.",
    timing: "4–8 weeks before",
  },
  {
    id: "passport-tickets",
    title: "Passport-ready attraction tickets",
    detail:
      "Forbidden City, some museums, and popular parks require real-name booking with your passport number. Book as soon as your dates firm up.",
    timing: "1–3 weeks before",
  },
  {
    id: "payments",
    title: "Alipay / WeChat Pay setup",
    detail:
      "Link an eligible foreign card in Alipay and WeChat before you fly. Test a tiny payment. Carry a small cash buffer for edge cases.",
    timing: "1–2 weeks before",
  },
  {
    id: "connectivity",
    title: "eSIM or travel SIM",
    detail:
      "Buy a China-capable eSIM (or plan a SIM on arrival). Do not rely on hotel Wi‑Fi alone for maps, DiDi, and tickets.",
    timing: "Before departure",
  },
  {
    id: "apps",
    title: "Install the survival apps",
    detail:
      "Alipay, WeChat, DiDi, a translation app, Trip.com or 12306 for trains, and Apple Maps / Amap. Log in while you still have reliable home internet.",
    timing: "Before departure",
  },
  {
    id: "rail",
    title: "Book long-distance trains",
    detail:
      "High-speed seats sell out around holidays. Prefer Trip.com in English if 12306 feels opaque — same trains, clearer UI.",
    timing: "As soon as dates lock",
  },
];
