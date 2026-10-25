import Image from "next/image";
import Link from "next/link";
import { destinations } from "@/content/destinations";

export default function HomePage() {
  return (
    <>
      <section className="hero">
        <div className="hero__media" aria-hidden>
          <Image
            src="/images/great-wall.jpg"
            alt=""
            fill
            priority
            sizes="100vw"
          />
          <div className="hero__shade" />
        </div>
        <div className="wrap">
          <div className="hero__content">
            <p className="eyebrow">China · for first-time visitors</p>
            <h1>EASTBOUND</h1>
            <p className="lede">
              Get in, pay, move, and follow a clear plan — English guides for
              the parts of China travel that usually go wrong.
            </p>
            <div className="hero__actions">
              <Link href="/plan" className="btn btn--primary">
                Start planning
              </Link>
              <Link href="/itineraries/classic-10" className="btn btn--ghost">
                See a 10-day route
              </Link>
            </div>
          </div>
        </div>
      </section>

      <section className="section">
        <div className="wrap">
          <div className="section-head">
            <p className="eyebrow">Where to begin</p>
            <h2>Four cities that teach you the country</h2>
            <p className="muted">
              Start with hubs that have good trains, clear sights, and enough
              English support to build confidence.
            </p>
          </div>
          <div className="destination-rail">
            {destinations.map((d) => (
              <Link key={d.slug} href={`/destinations/${d.slug}`} className="dest-link">
                <Image
                  src={d.image}
                  alt={d.imageAlt}
                  fill
                  sizes="(max-width: 720px) 100vw, 25vw"
                />
                <div className="dest-link__body">
                  <p className="dest-link__meta">
                    {d.region} · {d.days}
                  </p>
                  <h3>{d.name}</h3>
                  <p className="dest-link__meta">{d.bestFor}</p>
                </div>
              </Link>
            ))}
          </div>
        </div>
      </section>

      <section className="section section--tight">
        <div className="wrap split">
          <div>
            <p className="eyebrow">Before you fly</p>
            <h2 className="prose-title" style={{ fontSize: "clamp(1.8rem, 3vw, 2.6rem)" }}>
              Clear the gates that stop most trips
            </h2>
            <p className="muted">
              Visa rules, foreigner payment setup, eSIM, and attraction booking
              with your passport — sorted in one checklist.
            </p>
            <p style={{ marginTop: "1.5rem" }}>
              <Link href="/plan" className="btn btn--primary">
                Open the plan checklist
              </Link>
            </p>
          </div>
          <div className="panel stack-gaps">
            <div>
              <h3>Pay without panic</h3>
              <p className="muted">
                Alipay and WeChat Pay step-by-steps for foreign cards — tested
                before you land.
              </p>
            </div>
            <div>
              <h3>Ride the rail</h3>
              <p className="muted">
                High-speed booking in English, metro QR codes, and when DiDi
                beats a taxi queue.
              </p>
            </div>
            <div>
              <h3>Book the hard tickets</h3>
              <p className="muted">
                Passport-linked museums and Wall routes with timing that
                survives jet lag.
              </p>
            </div>
            <p>
              <Link href="/essentials">Essentials guide →</Link>
            </p>
          </div>
        </div>
      </section>

      <section className="section">
        <div className="wrap">
          <div className="cta-band">
            <h2>A classic first route, already mapped</h2>
            <p>
              Beijing → Xi&apos;an → Shanghai in ten days: history, Wall,
              Warriors, then a soft landing.
            </p>
            <Link href="/itineraries/classic-10" className="btn btn--primary">
              View the itinerary
            </Link>
          </div>
        </div>
      </section>
    </>
  );
}
